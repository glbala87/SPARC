import { useState, useEffect, useRef } from 'react'
import { Activity, TrendingUp, AlertTriangle } from 'lucide-react'

interface QCData {
  total_reads: number
  valid_barcodes: number
  valid_barcode_pct: number
  corrected_barcodes: number
  mapped_reads: number
  assigned_reads: number
  cells: number
  genes: number
  median_genes_per_cell: number
  median_umi_per_cell: number
  timestamp: number
}

interface QCDashboardProps {
  jobId: string
}

export function QCDashboard({ jobId }: QCDashboardProps) {
  const [qcData, setQcData] = useState<QCData | null>(null)
  const [history, setHistory] = useState<QCData[]>([])
  const wsRef = useRef<WebSocket | null>(null)

  useEffect(() => {
    if (!jobId) return

    const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:'
    const wsUrl = `${protocol}//${window.location.host}/ws/pipeline/${jobId}`

    const ws = new WebSocket(wsUrl)
    wsRef.current = ws

    ws.onmessage = (event) => {
      try {
        const data = JSON.parse(event.data)
        if (data.type === 'qc_update' && data.qc_data) {
          const qc = { ...data.qc_data, timestamp: Date.now() }
          setQcData(qc)
          setHistory((prev) => [...prev.slice(-50), qc])
        }
      } catch (e) {
        console.error('Failed to parse QC update:', e)
      }
    }

    return () => {
      ws.close()
    }
  }, [jobId])

  if (!qcData) {
    return (
      <div className="bg-white shadow rounded-lg p-6">
        <div className="flex items-center space-x-2 text-gray-500">
          <Activity className="h-5 w-5 animate-pulse" />
          <span>Waiting for QC data...</span>
        </div>
      </div>
    )
  }

  const metrics = [
    { label: 'Total Reads', value: qcData.total_reads?.toLocaleString() || '0' },
    { label: 'Valid Barcodes', value: `${qcData.valid_barcodes?.toLocaleString() || '0'} (${qcData.valid_barcode_pct?.toFixed(1) || '0'}%)` },
    { label: 'Mapped Reads', value: qcData.mapped_reads?.toLocaleString() || '0' },
    { label: 'Assigned Reads', value: qcData.assigned_reads?.toLocaleString() || '0' },
    { label: 'Cells Detected', value: qcData.cells?.toLocaleString() || '0' },
    { label: 'Genes Detected', value: qcData.genes?.toLocaleString() || '0' },
    { label: 'Median Genes/Cell', value: qcData.median_genes_per_cell?.toFixed(0) || '0' },
    { label: 'Median UMIs/Cell', value: qcData.median_umi_per_cell?.toFixed(0) || '0' },
  ]

  const warnings: string[] = []
  if (qcData.valid_barcode_pct < 50) warnings.push('Low barcode match rate')
  if (qcData.median_genes_per_cell < 500) warnings.push('Low gene detection per cell')
  if (qcData.cells < 100) warnings.push('Low cell count')

  return (
    <div className="space-y-4">
      <div className="bg-white shadow rounded-lg p-6">
        <div className="flex items-center justify-between mb-4">
          <h3 className="text-lg font-semibold text-gray-900 flex items-center space-x-2">
            <TrendingUp className="h-5 w-5 text-indigo-500" />
            <span>Real-time QC Dashboard</span>
          </h3>
          <span className="text-xs text-gray-400">
            Updates: {history.length}
          </span>
        </div>

        <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
          {metrics.map((m) => (
            <div key={m.label} className="bg-gray-50 rounded-lg p-3">
              <dt className="text-xs text-gray-500 uppercase tracking-wide">{m.label}</dt>
              <dd className="text-lg font-semibold text-gray-900 mt-1">{m.value}</dd>
            </div>
          ))}
        </div>
      </div>

      {warnings.length > 0 && (
        <div className="bg-yellow-50 border border-yellow-200 rounded-lg p-4">
          <div className="flex items-center space-x-2 text-yellow-800 mb-2">
            <AlertTriangle className="h-5 w-5" />
            <span className="font-medium">QC Warnings</span>
          </div>
          <ul className="list-disc list-inside text-sm text-yellow-700 space-y-1">
            {warnings.map((w, i) => (
              <li key={i}>{w}</li>
            ))}
          </ul>
        </div>
      )}
    </div>
  )
}
