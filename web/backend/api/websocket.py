"""
WebSocket handlers for real-time progress updates.
"""

import asyncio
from typing import Dict, Set

from fastapi import APIRouter, WebSocket, WebSocketDisconnect

websocket_router = APIRouter()

# Connected clients per job
connected_clients: Dict[str, Set[WebSocket]] = {}


class ConnectionManager:
    """Manage WebSocket connections."""

    def __init__(self):
        self.active_connections: Dict[str, Set[WebSocket]] = {}

    async def connect(self, websocket: WebSocket, job_id: str):
        """Accept and register a WebSocket connection."""
        await websocket.accept()
        if job_id not in self.active_connections:
            self.active_connections[job_id] = set()
        self.active_connections[job_id].add(websocket)

    def disconnect(self, websocket: WebSocket, job_id: str):
        """Remove a WebSocket connection."""
        if job_id in self.active_connections:
            self.active_connections[job_id].discard(websocket)
            if not self.active_connections[job_id]:
                del self.active_connections[job_id]

    async def broadcast(self, job_id: str, message: dict):
        """Broadcast a message to all clients for a job."""
        if job_id in self.active_connections:
            dead_connections = set()
            for connection in self.active_connections[job_id]:
                try:
                    await connection.send_json(message)
                except Exception:
                    dead_connections.add(connection)

            # Clean up dead connections
            for connection in dead_connections:
                self.active_connections[job_id].discard(connection)


manager = ConnectionManager()


@websocket_router.websocket("/pipeline/{job_id}")
async def websocket_pipeline_progress(websocket: WebSocket, job_id: str):
    """WebSocket endpoint for pipeline progress updates."""
    await manager.connect(websocket, job_id)

    try:
        while True:
            # Wait for messages from client (heartbeat)
            try:
                data = await asyncio.wait_for(websocket.receive_json(), timeout=30.0)

                if data.get("type") == "ping":
                    await websocket.send_json({"type": "pong"})

            except asyncio.TimeoutError:
                # Send heartbeat
                await websocket.send_json({"type": "heartbeat"})

    except WebSocketDisconnect:
        manager.disconnect(websocket, job_id)


async def send_progress_update(job_id: str, progress: float, message: str, status: str = "running"):
    """Send a progress update to all connected clients."""
    await manager.broadcast(
        job_id,
        {
            "type": "progress",
            "job_id": job_id,
            "progress": progress,
            "message": message,
            "status": status,
        },
    )


async def send_result(job_id: str, result: dict):
    """Send pipeline result to all connected clients."""
    await manager.broadcast(
        job_id,
        {
            "type": "result",
            "job_id": job_id,
            "result": result,
        },
    )
