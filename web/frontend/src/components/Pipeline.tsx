import { useState } from 'react'
import { useMutation, useQuery } from '@tanstack/react-query'
import { Play, Loader2, CheckCircle, XCircle } from 'lucide-react'
import axios from 'axios'
import { useWebSocket } from '../hooks/useWebSocket'

interface Protocol {
  id: string
  name: string
  barcode_len: number
  umi_len: number
}

interface PipelineConfig {
  sample_name: string
  protocol: string
  max_mismatch: number
  min_genes: number
  max_genes: number
  max_mito: number
  resolution: number
}

export function Pipeline() {
  const [jobId, setJobId] = useState('')
  const [config, setConfig] = useState<PipelineConfig>({
    sample_name: 'sample',
    protocol: '10x-3prime-v3',
    max_mismatch: 1,
    min_genes: 200,
    max_genes: 10000,
    max_mito: 20.0,
    resolution: 1.0,
  })

  const { data: protocols } = useQuery({
    queryKey: ['protocols'],
    queryFn: async () => {
      const response = await axios.get<{ protocols: Protocol[] }>('/api/protocols')
      return response.data.protocols
    },
  })

  const { status, progress, message } = useWebSocket(jobId)

  const startPipeline = useMutation({
    mutationFn: async () => {
      const response = await axios.post(`/api/pipeline/${jobId}`, config)
      return response.data
    },
  })

  const statusQuery = useQuery({
    queryKey: ['pipeline-status', jobId],
    queryFn: async () => {
      const response = await axios.get(`/api/pipeline/${jobId}/status`)
      return response.data
    },
    enabled: !!jobId && startPipeline.isSuccess,
    refetchInterval: status === 'running' ? 2000 : false,
  })

  return (
    <div className="space-y-6">
      <div className="bg-white shadow rounded-lg p-6">
        <h2 className="text-2xl font-bold text-gray-900 mb-4">Pipeline Configuration</h2>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          {/* Job ID */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Job ID
            </label>
            <input
              type="text"
              value={jobId}
              onChange={(e) => setJobId(e.target.value)}
              placeholder="Enter job ID from upload"
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-indigo-500"
            />
          </div>

          {/* Sample Name */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Sample Name
            </label>
            <input
              type="text"
              value={config.sample_name}
              onChange={(e) => setConfig({ ...config, sample_name: e.target.value })}
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-indigo-500"
            />
          </div>

          {/* Protocol */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Protocol
            </label>
            <select
              value={config.protocol}
              onChange={(e) => setConfig({ ...config, protocol: e.target.value })}
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-indigo-500"
            >
              {protocols?.map((p) => (
                <option key={p.id} value={p.id}>
                  {p.name} (BC: {p.barcode_len}bp, UMI: {p.umi_len}bp)
                </option>
              ))}
            </select>
          </div>

          {/* Max Mismatch */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Max Barcode Mismatch
            </label>
            <input
              type="number"
              min={0}
              max={3}
              value={config.max_mismatch}
              onChange={(e) => setConfig({ ...config, max_mismatch: parseInt(e.target.value) })}
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-indigo-500"
            />
          </div>

          {/* Min Genes */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Min Genes per Cell
            </label>
            <input
              type="number"
              min={0}
              value={config.min_genes}
              onChange={(e) => setConfig({ ...config, min_genes: parseInt(e.target.value) })}
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-indigo-500"
            />
          </div>

          {/* Max Genes */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Max Genes per Cell
            </label>
            <input
              type="number"
              min={0}
              value={config.max_genes}
              onChange={(e) => setConfig({ ...config, max_genes: parseInt(e.target.value) })}
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-indigo-500"
            />
          </div>

          {/* Max Mito */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Max Mitochondrial %
            </label>
            <input
              type="number"
              min={0}
              max={100}
              step={0.1}
              value={config.max_mito}
              onChange={(e) => setConfig({ ...config, max_mito: parseFloat(e.target.value) })}
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-indigo-500"
            />
          </div>

          {/* Resolution */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Clustering Resolution
            </label>
            <input
              type="number"
              min={0.1}
              max={3}
              step={0.1}
              value={config.resolution}
              onChange={(e) => setConfig({ ...config, resolution: parseFloat(e.target.value) })}
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-indigo-500"
            />
          </div>
        </div>

        <div className="mt-6 flex justify-end">
          <button
            onClick={() => startPipeline.mutate()}
            disabled={!jobId || startPipeline.isPending}
            className="bg-indigo-600 text-white px-6 py-2 rounded-md hover:bg-indigo-700
              disabled:bg-gray-400 disabled:cursor-not-allowed flex items-center space-x-2"
          >
            {startPipeline.isPending ? (
              <>
                <Loader2 className="h-4 w-4 animate-spin" />
                <span>Starting...</span>
              </>
            ) : (
              <>
                <Play className="h-4 w-4" />
                <span>Start Pipeline</span>
              </>
            )}
          </button>
        </div>
      </div>

      {/* Progress */}
      {(startPipeline.isSuccess || statusQuery.data) && (
        <div className="bg-white shadow rounded-lg p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">Pipeline Progress</h3>

          <div className="space-y-4">
            <div className="flex items-center space-x-4">
              {status === 'completed' ? (
                <CheckCircle className="h-6 w-6 text-green-500" />
              ) : status === 'failed' ? (
                <XCircle className="h-6 w-6 text-red-500" />
              ) : (
                <Loader2 className="h-6 w-6 text-indigo-500 animate-spin" />
              )}
              <span className="text-gray-700">{message || statusQuery.data?.message}</span>
            </div>

            <div className="w-full bg-gray-200 rounded-full h-2.5">
              <div
                className="bg-indigo-600 h-2.5 rounded-full transition-all duration-500"
                style={{ width: `${(progress || statusQuery.data?.progress || 0) * 100}%` }}
              />
            </div>

            <p className="text-sm text-gray-500">
              {Math.round((progress || statusQuery.data?.progress || 0) * 100)}% complete
            </p>
          </div>

          {statusQuery.data?.result && (
            <div className="mt-6 p-4 bg-gray-50 rounded-md">
              <h4 className="font-medium text-gray-900 mb-2">Results</h4>
              <dl className="grid grid-cols-2 gap-4 text-sm">
                <div>
                  <dt className="text-gray-500">Total Reads</dt>
                  <dd className="font-medium">{statusQuery.data.result.total_reads?.toLocaleString()}</dd>
                </div>
                <div>
                  <dt className="text-gray-500">Valid Barcodes</dt>
                  <dd className="font-medium">{statusQuery.data.result.valid_barcodes?.toLocaleString()}</dd>
                </div>
                <div>
                  <dt className="text-gray-500">Cells</dt>
                  <dd className="font-medium">{statusQuery.data.result.cells?.toLocaleString()}</dd>
                </div>
                <div>
                  <dt className="text-gray-500">Genes</dt>
                  <dd className="font-medium">{statusQuery.data.result.genes?.toLocaleString()}</dd>
                </div>
              </dl>
            </div>
          )}
        </div>
      )}
    </div>
  )
}
