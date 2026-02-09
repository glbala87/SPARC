import { useEffect, useState, useRef, useCallback } from 'react'

interface WebSocketMessage {
  type: 'progress' | 'result' | 'heartbeat' | 'pong'
  job_id?: string
  progress?: number
  message?: string
  status?: string
  result?: Record<string, unknown>
}

interface UseWebSocketReturn {
  status: string
  progress: number
  message: string
  result: Record<string, unknown> | null
  isConnected: boolean
}

export function useWebSocket(jobId: string | null): UseWebSocketReturn {
  const [status, setStatus] = useState('pending')
  const [progress, setProgress] = useState(0)
  const [message, setMessage] = useState('')
  const [result, setResult] = useState<Record<string, unknown> | null>(null)
  const [isConnected, setIsConnected] = useState(false)

  const wsRef = useRef<WebSocket | null>(null)
  const reconnectTimeoutRef = useRef<number | null>(null)

  const connect = useCallback(() => {
    if (!jobId) return

    const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:'
    const wsUrl = `${protocol}//${window.location.host}/ws/pipeline/${jobId}`

    try {
      const ws = new WebSocket(wsUrl)
      wsRef.current = ws

      ws.onopen = () => {
        setIsConnected(true)
        console.log('WebSocket connected')
      }

      ws.onmessage = (event) => {
        try {
          const data: WebSocketMessage = JSON.parse(event.data)

          switch (data.type) {
            case 'progress':
              if (data.progress !== undefined) setProgress(data.progress)
              if (data.message) setMessage(data.message)
              if (data.status) setStatus(data.status)
              break
            case 'result':
              if (data.result) setResult(data.result)
              setStatus('completed')
              break
            case 'heartbeat':
              // Respond to heartbeat
              ws.send(JSON.stringify({ type: 'ping' }))
              break
            case 'pong':
              // Heartbeat acknowledged
              break
          }
        } catch (e) {
          console.error('Failed to parse WebSocket message:', e)
        }
      }

      ws.onerror = (error) => {
        console.error('WebSocket error:', error)
      }

      ws.onclose = () => {
        setIsConnected(false)
        console.log('WebSocket disconnected')

        // Reconnect after 5 seconds if job is still in progress
        if (status === 'running') {
          reconnectTimeoutRef.current = window.setTimeout(() => {
            connect()
          }, 5000)
        }
      }
    } catch (e) {
      console.error('Failed to create WebSocket:', e)
    }
  }, [jobId, status])

  useEffect(() => {
    connect()

    return () => {
      if (wsRef.current) {
        wsRef.current.close()
      }
      if (reconnectTimeoutRef.current) {
        clearTimeout(reconnectTimeoutRef.current)
      }
    }
  }, [connect])

  return { status, progress, message, result, isConnected }
}
