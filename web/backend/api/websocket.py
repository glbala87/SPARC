"""
WebSocket handlers for real-time progress updates.
"""

import asyncio
import logging
import os
from typing import Dict, Set

from fastapi import APIRouter, WebSocket, WebSocketDisconnect

logger = logging.getLogger("sparc.websocket")

websocket_router = APIRouter()

MAX_CONNECTIONS_PER_JOB = int(os.getenv("SPARC_WS_MAX_PER_JOB", "10"))
MAX_TOTAL_CONNECTIONS = int(os.getenv("SPARC_WS_MAX_TOTAL", "100"))


class ConnectionManager:
    """Manage WebSocket connections with limits and thread safety."""

    def __init__(self):
        self.active_connections: Dict[str, Set[WebSocket]] = {}
        self._lock = asyncio.Lock()

    @property
    def total_connections(self) -> int:
        return sum(len(conns) for conns in self.active_connections.values())

    async def connect(self, websocket: WebSocket, job_id: str) -> bool:
        """Accept and register a WebSocket connection. Returns False if limit exceeded."""
        async with self._lock:
            if self.total_connections >= MAX_TOTAL_CONNECTIONS:
                logger.warning("WebSocket connection rejected: total limit reached (%d)", MAX_TOTAL_CONNECTIONS)
                await websocket.close(code=1013, reason="Server at capacity")
                return False

            job_conns = self.active_connections.get(job_id, set())
            if len(job_conns) >= MAX_CONNECTIONS_PER_JOB:
                logger.warning("WebSocket connection rejected for job %s: per-job limit reached (%d)",
                               job_id, MAX_CONNECTIONS_PER_JOB)
                await websocket.close(code=1013, reason="Too many connections for this job")
                return False

            await websocket.accept()
            if job_id not in self.active_connections:
                self.active_connections[job_id] = set()
            self.active_connections[job_id].add(websocket)

        logger.info("WebSocket connected for job %s (total=%d)", job_id, self.total_connections)
        return True

    async def disconnect(self, websocket: WebSocket, job_id: str):
        """Remove a WebSocket connection."""
        async with self._lock:
            if job_id in self.active_connections:
                self.active_connections[job_id].discard(websocket)
                if not self.active_connections[job_id]:
                    del self.active_connections[job_id]
        logger.info("WebSocket disconnected for job %s (total=%d)", job_id, self.total_connections)

    async def broadcast(self, job_id: str, message: dict):
        """Broadcast a message to all clients for a job."""
        async with self._lock:
            connections = self.active_connections.get(job_id, set()).copy()

        dead_connections = set()
        for connection in connections:
            try:
                await connection.send_json(message)
            except Exception:
                dead_connections.add(connection)

        if dead_connections:
            async with self._lock:
                if job_id in self.active_connections:
                    self.active_connections[job_id] -= dead_connections
                    if not self.active_connections[job_id]:
                        del self.active_connections[job_id]

    async def cleanup_all(self):
        """Close all connections (for graceful shutdown)."""
        async with self._lock:
            all_conns = [
                (job_id, ws)
                for job_id, conns in self.active_connections.items()
                for ws in conns
            ]
        for job_id, ws in all_conns:
            try:
                await ws.close(code=1001, reason="Server shutting down")
            except Exception:
                pass
        async with self._lock:
            self.active_connections.clear()
        logger.info("All WebSocket connections closed")


manager = ConnectionManager()


@websocket_router.websocket("/pipeline/{job_id}")
async def websocket_pipeline_progress(websocket: WebSocket, job_id: str):
    """WebSocket endpoint for pipeline progress updates."""
    connected = await manager.connect(websocket, job_id)
    if not connected:
        return

    try:
        while True:
            try:
                data = await asyncio.wait_for(websocket.receive_json(), timeout=30.0)
                if data.get("type") == "ping":
                    await websocket.send_json({"type": "pong"})
            except asyncio.TimeoutError:
                await websocket.send_json({"type": "heartbeat"})
    except WebSocketDisconnect:
        await manager.disconnect(websocket, job_id)
    except Exception:
        await manager.disconnect(websocket, job_id)


async def send_progress_update(job_id: str, progress: float, message: str, status: str = "running"):
    """Send a progress update to all connected clients."""
    await manager.broadcast(job_id, {
        "type": "progress", "job_id": job_id,
        "progress": progress, "message": message, "status": status,
    })


async def send_result(job_id: str, result: dict):
    """Send pipeline result to all connected clients."""
    await manager.broadcast(job_id, {"type": "result", "job_id": job_id, "result": result})


async def send_qc_update(job_id: str, qc_data: dict):
    """Send real-time QC metrics update to all connected clients."""
    await manager.broadcast(job_id, {"type": "qc_update", "job_id": job_id, "qc_data": qc_data})
