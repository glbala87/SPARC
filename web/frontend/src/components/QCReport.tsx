import { useState } from 'react'
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  ScatterChart,
  Scatter,
} from 'recharts'
import { AlertTriangle, CheckCircle } from 'lucide-react'

// Mock data for demonstration
const mockQCData = {
  sample_name: 'Sample1',
  metrics: {
    total_reads: 10000000,
    valid_barcode_reads: 9500000,
    mapped_reads: 9000000,
    assigned_reads: 8500000,
    num_cells: 5000,
    mean_reads_per_cell: 1900,
    median_reads_per_cell: 1800,
    mean_genes_per_cell: 2500,
    median_genes_per_cell: 2400,
    total_genes: 20000,
    mean_umi_per_cell: 5000,
    median_umi_per_cell: 4800,
    sequencing_saturation: 0.65,
    fraction_reads_in_cells: 0.85,
  },
  warnings: ['Low sequencing saturation (65%)'],
  cells_per_gene_histogram: [
    { range: '1-10', count: 5000 },
    { range: '11-50', count: 4000 },
    { range: '51-100', count: 3000 },
    { range: '101-500', count: 5000 },
    { range: '501-1000', count: 2000 },
    { range: '>1000', count: 1000 },
  ],
  genes_per_cell_histogram: [
    { range: '0-500', count: 200 },
    { range: '500-1000', count: 500 },
    { range: '1000-2000', count: 1500 },
    { range: '2000-3000', count: 1800 },
    { range: '3000-4000', count: 700 },
    { range: '>4000', count: 300 },
  ],
  scatter_data: Array.from({ length: 100 }, (_, i) => ({
    total_counts: Math.random() * 10000 + 1000,
    n_genes: Math.random() * 3000 + 500,
    mito_percent: Math.random() * 20,
  })),
}

export function QCReport() {
  const [data] = useState(mockQCData)

  return (
    <div className="space-y-6">
      {/* Summary Stats */}
      <div className="bg-white shadow rounded-lg p-6">
        <h2 className="text-2xl font-bold text-gray-900 mb-4">QC Summary</h2>

        <div className="grid grid-cols-2 md:grid-cols-4 gap-6">
          <StatCard
            label="Total Cells"
            value={data.metrics.num_cells.toLocaleString()}
          />
          <StatCard
            label="Median Genes/Cell"
            value={data.metrics.median_genes_per_cell.toLocaleString()}
          />
          <StatCard
            label="Median UMIs/Cell"
            value={data.metrics.median_umi_per_cell.toLocaleString()}
          />
          <StatCard
            label="Sequencing Saturation"
            value={`${(data.metrics.sequencing_saturation * 100).toFixed(1)}%`}
          />
          <StatCard
            label="Total Reads"
            value={`${(data.metrics.total_reads / 1e6).toFixed(1)}M`}
          />
          <StatCard
            label="Valid Barcodes"
            value={`${((data.metrics.valid_barcode_reads / data.metrics.total_reads) * 100).toFixed(1)}%`}
          />
          <StatCard
            label="Mapping Rate"
            value={`${((data.metrics.mapped_reads / data.metrics.total_reads) * 100).toFixed(1)}%`}
          />
          <StatCard
            label="Reads in Cells"
            value={`${(data.metrics.fraction_reads_in_cells * 100).toFixed(1)}%`}
          />
        </div>

        {/* Warnings */}
        {data.warnings.length > 0 && (
          <div className="mt-6 p-4 bg-yellow-50 rounded-md">
            <div className="flex items-center space-x-2 text-yellow-700">
              <AlertTriangle className="h-5 w-5" />
              <span className="font-medium">Warnings</span>
            </div>
            <ul className="mt-2 list-disc list-inside text-yellow-600">
              {data.warnings.map((warning, i) => (
                <li key={i}>{warning}</li>
              ))}
            </ul>
          </div>
        )}

        {data.warnings.length === 0 && (
          <div className="mt-6 p-4 bg-green-50 rounded-md">
            <div className="flex items-center space-x-2 text-green-700">
              <CheckCircle className="h-5 w-5" />
              <span className="font-medium">All QC metrics passed</span>
            </div>
          </div>
        )}
      </div>

      {/* Charts */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        {/* Genes per Cell Histogram */}
        <div className="bg-white shadow rounded-lg p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">Genes per Cell</h3>
          <ResponsiveContainer width="100%" height={300}>
            <BarChart data={data.genes_per_cell_histogram}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis dataKey="range" />
              <YAxis />
              <Tooltip />
              <Bar dataKey="count" fill="#6366f1" />
            </BarChart>
          </ResponsiveContainer>
        </div>

        {/* Cells per Gene Histogram */}
        <div className="bg-white shadow rounded-lg p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">Cells per Gene</h3>
          <ResponsiveContainer width="100%" height={300}>
            <BarChart data={data.cells_per_gene_histogram}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis dataKey="range" />
              <YAxis />
              <Tooltip />
              <Bar dataKey="count" fill="#10b981" />
            </BarChart>
          </ResponsiveContainer>
        </div>

        {/* Scatter Plot */}
        <div className="bg-white shadow rounded-lg p-6 md:col-span-2">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">
            Total Counts vs Genes Detected
          </h3>
          <ResponsiveContainer width="100%" height={400}>
            <ScatterChart>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis type="number" dataKey="total_counts" name="Total Counts" />
              <YAxis type="number" dataKey="n_genes" name="Genes" />
              <Tooltip cursor={{ strokeDasharray: '3 3' }} />
              <Scatter
                name="Cells"
                data={data.scatter_data}
                fill="#6366f1"
                fillOpacity={0.6}
              />
            </ScatterChart>
          </ResponsiveContainer>
        </div>
      </div>
    </div>
  )
}

function StatCard({ label, value }: { label: string; value: string }) {
  return (
    <div className="bg-gray-50 rounded-lg p-4">
      <p className="text-sm text-gray-500">{label}</p>
      <p className="text-2xl font-bold text-gray-900">{value}</p>
    </div>
  )
}
