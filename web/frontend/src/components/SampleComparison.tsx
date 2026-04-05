import { useState } from 'react'
import { useQuery } from '@tanstack/react-query'
import { BarChart2, Plus, X, Loader2 } from 'lucide-react'
import axios from 'axios'

interface SampleResult {
  job_id: string
  sample_name: string
  total_reads: number
  valid_barcodes: number
  cells: number
  genes: number
  median_genes_per_cell: number
}

export function SampleComparison() {
  const [jobIds, setJobIds] = useState<string[]>([''])

  const addJobId = () => setJobIds([...jobIds, ''])
  const removeJobId = (index: number) => setJobIds(jobIds.filter((_, i) => i !== index))
  const updateJobId = (index: number, value: string) => {
    const updated = [...jobIds]
    updated[index] = value
    setJobIds(updated)
  }

  const validIds = jobIds.filter((id) => id.trim().length > 0)

  const { data, isLoading, refetch } = useQuery({
    queryKey: ['compare', validIds],
    queryFn: async () => {
      const response = await axios.post<{ samples: SampleResult[] }>('/api/compare', validIds)
      return response.data.samples
    },
    enabled: false,
  })

  return (
    <div className="space-y-6">
      <div className="bg-white shadow rounded-lg p-6">
        <h2 className="text-2xl font-bold text-gray-900 mb-4 flex items-center space-x-2">
          <BarChart2 className="h-6 w-6 text-indigo-500" />
          <span>Multi-Sample Comparison</span>
        </h2>

        <div className="space-y-3">
          {jobIds.map((id, index) => (
            <div key={index} className="flex items-center space-x-2">
              <input
                type="text"
                value={id}
                onChange={(e) => updateJobId(index, e.target.value)}
                placeholder={`Job ID #${index + 1}`}
                className="flex-1 px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-indigo-500"
              />
              {jobIds.length > 1 && (
                <button
                  onClick={() => removeJobId(index)}
                  className="p-2 text-red-500 hover:bg-red-50 rounded-md"
                >
                  <X className="h-4 w-4" />
                </button>
              )}
            </div>
          ))}

          <div className="flex space-x-3">
            <button
              onClick={addJobId}
              className="flex items-center space-x-1 px-3 py-2 border border-gray-300 rounded-md hover:bg-gray-50"
            >
              <Plus className="h-4 w-4" />
              <span>Add Sample</span>
            </button>

            <button
              onClick={() => refetch()}
              disabled={validIds.length < 2 || isLoading}
              className="bg-indigo-600 text-white px-4 py-2 rounded-md hover:bg-indigo-700
                disabled:bg-gray-400 disabled:cursor-not-allowed flex items-center space-x-2"
            >
              {isLoading ? (
                <Loader2 className="h-4 w-4 animate-spin" />
              ) : (
                <BarChart2 className="h-4 w-4" />
              )}
              <span>Compare</span>
            </button>
          </div>
        </div>
      </div>

      {data && data.length > 0 && (
        <div className="bg-white shadow rounded-lg p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">Comparison Results</h3>

          <div className="overflow-x-auto">
            <table className="min-w-full divide-y divide-gray-200">
              <thead className="bg-gray-50">
                <tr>
                  <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">Sample</th>
                  <th className="px-4 py-3 text-right text-xs font-medium text-gray-500 uppercase">Total Reads</th>
                  <th className="px-4 py-3 text-right text-xs font-medium text-gray-500 uppercase">Valid Barcodes</th>
                  <th className="px-4 py-3 text-right text-xs font-medium text-gray-500 uppercase">Cells</th>
                  <th className="px-4 py-3 text-right text-xs font-medium text-gray-500 uppercase">Genes</th>
                  <th className="px-4 py-3 text-right text-xs font-medium text-gray-500 uppercase">Median Genes/Cell</th>
                </tr>
              </thead>
              <tbody className="bg-white divide-y divide-gray-200">
                {data.map((sample) => (
                  <tr key={sample.job_id}>
                    <td className="px-4 py-3 text-sm font-medium text-gray-900">{sample.sample_name}</td>
                    <td className="px-4 py-3 text-sm text-gray-700 text-right">{sample.total_reads?.toLocaleString()}</td>
                    <td className="px-4 py-3 text-sm text-gray-700 text-right">{sample.valid_barcodes?.toLocaleString()}</td>
                    <td className="px-4 py-3 text-sm text-gray-700 text-right">{sample.cells?.toLocaleString()}</td>
                    <td className="px-4 py-3 text-sm text-gray-700 text-right">{sample.genes?.toLocaleString()}</td>
                    <td className="px-4 py-3 text-sm text-gray-700 text-right">{sample.median_genes_per_cell?.toLocaleString()}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}
    </div>
  )
}
