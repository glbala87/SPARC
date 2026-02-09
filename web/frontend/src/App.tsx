import { useState } from 'react'
import { BrowserRouter, Routes, Route, Link } from 'react-router-dom'
import { Upload } from './components/Upload'
import { Pipeline } from './components/Pipeline'
import { QCReport } from './components/QCReport'
import { Visualization } from './components/Visualization'
import { Dna, BarChart2, FileUp, Settings } from 'lucide-react'

function Navigation() {
  return (
    <nav className="bg-indigo-600 text-white shadow-lg">
      <div className="max-w-7xl mx-auto px-4">
        <div className="flex items-center justify-between h-16">
          <Link to="/" className="flex items-center space-x-2">
            <Dna className="h-8 w-8" />
            <span className="font-bold text-xl">SPARC</span>
          </Link>
          <div className="flex space-x-4">
            <Link
              to="/"
              className="flex items-center space-x-1 px-3 py-2 rounded-md hover:bg-indigo-500"
            >
              <FileUp className="h-4 w-4" />
              <span>Upload</span>
            </Link>
            <Link
              to="/pipeline"
              className="flex items-center space-x-1 px-3 py-2 rounded-md hover:bg-indigo-500"
            >
              <Settings className="h-4 w-4" />
              <span>Pipeline</span>
            </Link>
            <Link
              to="/qc"
              className="flex items-center space-x-1 px-3 py-2 rounded-md hover:bg-indigo-500"
            >
              <BarChart2 className="h-4 w-4" />
              <span>QC Report</span>
            </Link>
            <Link
              to="/visualization"
              className="flex items-center space-x-1 px-3 py-2 rounded-md hover:bg-indigo-500"
            >
              <Dna className="h-4 w-4" />
              <span>Visualization</span>
            </Link>
          </div>
        </div>
      </div>
    </nav>
  )
}

function App() {
  return (
    <BrowserRouter>
      <div className="min-h-screen bg-gray-50">
        <Navigation />
        <main className="max-w-7xl mx-auto py-6 px-4">
          <Routes>
            <Route path="/" element={<Upload />} />
            <Route path="/pipeline" element={<Pipeline />} />
            <Route path="/qc" element={<QCReport />} />
            <Route path="/visualization" element={<Visualization />} />
          </Routes>
        </main>
      </div>
    </BrowserRouter>
  )
}

export default App
