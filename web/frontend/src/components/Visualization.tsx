import { useState, useMemo } from 'react'
import { ScatterChart, Scatter, XAxis, YAxis, Tooltip, ResponsiveContainer } from 'recharts'

// Generate mock UMAP data
function generateMockUMAP(numCells: number, numClusters: number) {
  const clusters: { x: number; y: number; cluster: number; genes: number; umis: number }[] = []

  for (let c = 0; c < numClusters; c++) {
    const centerX = (Math.random() - 0.5) * 20
    const centerY = (Math.random() - 0.5) * 20
    const cellsInCluster = Math.floor(numCells / numClusters)

    for (let i = 0; i < cellsInCluster; i++) {
      clusters.push({
        x: centerX + (Math.random() - 0.5) * 4,
        y: centerY + (Math.random() - 0.5) * 4,
        cluster: c,
        genes: Math.floor(Math.random() * 3000 + 500),
        umis: Math.floor(Math.random() * 10000 + 1000),
      })
    }
  }

  return clusters
}

const CLUSTER_COLORS = [
  '#e11d48', '#f97316', '#eab308', '#22c55e', '#14b8a6',
  '#0ea5e9', '#6366f1', '#a855f7', '#ec4899', '#78716c',
]

export function Visualization() {
  const [numClusters] = useState(8)
  const [colorBy, setColorBy] = useState<'cluster' | 'genes' | 'umis'>('cluster')
  const [selectedCluster, setSelectedCluster] = useState<number | null>(null)

  const data = useMemo(() => generateMockUMAP(5000, numClusters), [numClusters])

  const filteredData = useMemo(() => {
    if (selectedCluster === null) return data
    return data.filter((d) => d.cluster === selectedCluster)
  }, [data, selectedCluster])

  const clusterCounts = useMemo(() => {
    const counts: Record<number, number> = {}
    data.forEach((d) => {
      counts[d.cluster] = (counts[d.cluster] || 0) + 1
    })
    return counts
  }, [data])

  return (
    <div className="space-y-6">
      <div className="bg-white shadow rounded-lg p-6">
        <div className="flex items-center justify-between mb-4">
          <h2 className="text-2xl font-bold text-gray-900">UMAP Visualization</h2>
          <div className="flex items-center space-x-4">
            <label className="text-sm text-gray-600">Color by:</label>
            <select
              value={colorBy}
              onChange={(e) => setColorBy(e.target.value as 'cluster' | 'genes' | 'umis')}
              className="px-3 py-1 border border-gray-300 rounded-md text-sm"
            >
              <option value="cluster">Cluster</option>
              <option value="genes">Genes</option>
              <option value="umis">UMIs</option>
            </select>
          </div>
        </div>

        <div className="flex">
          {/* UMAP Plot */}
          <div className="flex-1">
            <ResponsiveContainer width="100%" height={500}>
              <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
                <XAxis type="number" dataKey="x" name="UMAP1" tick={false} axisLine={false} />
                <YAxis type="number" dataKey="y" name="UMAP2" tick={false} axisLine={false} />
                <Tooltip
                  content={({ active, payload }) => {
                    if (active && payload && payload.length) {
                      const d = payload[0].payload
                      return (
                        <div className="bg-white p-2 shadow rounded border text-sm">
                          <p>Cluster: {d.cluster}</p>
                          <p>Genes: {d.genes}</p>
                          <p>UMIs: {d.umis}</p>
                        </div>
                      )
                    }
                    return null
                  }}
                />
                {colorBy === 'cluster' ? (
                  // Render each cluster as separate scatter
                  Array.from(new Set(filteredData.map((d) => d.cluster))).map((cluster) => (
                    <Scatter
                      key={cluster}
                      name={`Cluster ${cluster}`}
                      data={filteredData.filter((d) => d.cluster === cluster)}
                      fill={CLUSTER_COLORS[cluster % CLUSTER_COLORS.length]}
                      fillOpacity={0.7}
                    />
                  ))
                ) : (
                  <Scatter
                    name="Cells"
                    data={filteredData}
                    fill="#6366f1"
                    fillOpacity={0.7}
                  />
                )}
              </ScatterChart>
            </ResponsiveContainer>
          </div>

          {/* Legend */}
          <div className="w-48 pl-4">
            <h3 className="text-sm font-medium text-gray-700 mb-2">Clusters</h3>
            <div className="space-y-1">
              <button
                onClick={() => setSelectedCluster(null)}
                className={`w-full text-left px-2 py-1 rounded text-sm ${
                  selectedCluster === null ? 'bg-gray-200' : 'hover:bg-gray-100'
                }`}
              >
                All ({data.length})
              </button>
              {Object.entries(clusterCounts).map(([cluster, count]) => (
                <button
                  key={cluster}
                  onClick={() =>
                    setSelectedCluster(
                      selectedCluster === parseInt(cluster) ? null : parseInt(cluster)
                    )
                  }
                  className={`w-full text-left px-2 py-1 rounded text-sm flex items-center space-x-2 ${
                    selectedCluster === parseInt(cluster) ? 'bg-gray-200' : 'hover:bg-gray-100'
                  }`}
                >
                  <span
                    className="w-3 h-3 rounded-full"
                    style={{
                      backgroundColor: CLUSTER_COLORS[parseInt(cluster) % CLUSTER_COLORS.length],
                    }}
                  />
                  <span>
                    Cluster {cluster} ({count})
                  </span>
                </button>
              ))}
            </div>
          </div>
        </div>
      </div>

      {/* Cluster Stats */}
      <div className="bg-white shadow rounded-lg p-6">
        <h3 className="text-lg font-semibold text-gray-900 mb-4">Cluster Statistics</h3>
        <div className="overflow-x-auto">
          <table className="min-w-full divide-y divide-gray-200">
            <thead className="bg-gray-50">
              <tr>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">
                  Cluster
                </th>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">
                  Cells
                </th>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">
                  % of Total
                </th>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">
                  Median Genes
                </th>
                <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">
                  Median UMIs
                </th>
              </tr>
            </thead>
            <tbody className="bg-white divide-y divide-gray-200">
              {Object.entries(clusterCounts).map(([cluster, count]) => {
                const clusterData = data.filter((d) => d.cluster === parseInt(cluster))
                const genes = clusterData.map((d) => d.genes).sort((a, b) => a - b)
                const umis = clusterData.map((d) => d.umis).sort((a, b) => a - b)
                const medianGenes = genes[Math.floor(genes.length / 2)]
                const medianUmis = umis[Math.floor(umis.length / 2)]

                return (
                  <tr key={cluster} className="hover:bg-gray-50">
                    <td className="px-6 py-4 whitespace-nowrap">
                      <div className="flex items-center space-x-2">
                        <span
                          className="w-3 h-3 rounded-full"
                          style={{
                            backgroundColor:
                              CLUSTER_COLORS[parseInt(cluster) % CLUSTER_COLORS.length],
                          }}
                        />
                        <span>Cluster {cluster}</span>
                      </div>
                    </td>
                    <td className="px-6 py-4 whitespace-nowrap">{count}</td>
                    <td className="px-6 py-4 whitespace-nowrap">
                      {((count / data.length) * 100).toFixed(1)}%
                    </td>
                    <td className="px-6 py-4 whitespace-nowrap">{medianGenes}</td>
                    <td className="px-6 py-4 whitespace-nowrap">{medianUmis}</td>
                  </tr>
                )
              })}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  )
}
