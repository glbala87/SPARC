import { useState, useCallback } from 'react'
import { useMutation } from '@tanstack/react-query'
import { Upload as UploadIcon, File, Check, AlertCircle } from 'lucide-react'
import axios from 'axios'

interface UploadResponse {
  job_id: string
  files: {
    r1: string
    r2: string | null
    whitelist: string | null
  }
}

export function Upload() {
  const [r1File, setR1File] = useState<File | null>(null)
  const [r2File, setR2File] = useState<File | null>(null)
  const [whitelistFile, setWhitelistFile] = useState<File | null>(null)
  const [jobId, setJobId] = useState<string | null>(null)

  const uploadMutation = useMutation({
    mutationFn: async () => {
      const formData = new FormData()
      if (r1File) formData.append('r1', r1File)
      if (r2File) formData.append('r2', r2File)
      if (whitelistFile) formData.append('whitelist', whitelistFile)

      const response = await axios.post<UploadResponse>('/api/upload', formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
      })
      return response.data
    },
    onSuccess: (data) => {
      setJobId(data.job_id)
    },
  })

  const handleDrop = useCallback(
    (e: React.DragEvent, setter: (file: File) => void) => {
      e.preventDefault()
      const file = e.dataTransfer.files[0]
      if (file) setter(file)
    },
    []
  )

  const handleDragOver = (e: React.DragEvent) => {
    e.preventDefault()
  }

  return (
    <div className="space-y-6">
      <div className="bg-white shadow rounded-lg p-6">
        <h2 className="text-2xl font-bold text-gray-900 mb-4">Upload Files</h2>
        <p className="text-gray-600 mb-6">
          Upload your FASTQ files and barcode whitelist to start the analysis pipeline.
        </p>

        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          {/* R1 Upload */}
          <div
            className={`border-2 border-dashed rounded-lg p-6 text-center cursor-pointer
              ${r1File ? 'border-green-400 bg-green-50' : 'border-gray-300 hover:border-indigo-400'}`}
            onDrop={(e) => handleDrop(e, setR1File)}
            onDragOver={handleDragOver}
            onClick={() => document.getElementById('r1-input')?.click()}
          >
            <input
              id="r1-input"
              type="file"
              className="hidden"
              accept=".fastq,.fq,.gz"
              onChange={(e) => e.target.files?.[0] && setR1File(e.target.files[0])}
            />
            {r1File ? (
              <>
                <Check className="mx-auto h-12 w-12 text-green-500" />
                <p className="mt-2 text-sm text-green-600">{r1File.name}</p>
              </>
            ) : (
              <>
                <UploadIcon className="mx-auto h-12 w-12 text-gray-400" />
                <p className="mt-2 text-sm text-gray-600">R1 FASTQ (required)</p>
                <p className="text-xs text-gray-500">Drag & drop or click to upload</p>
              </>
            )}
          </div>

          {/* R2 Upload */}
          <div
            className={`border-2 border-dashed rounded-lg p-6 text-center cursor-pointer
              ${r2File ? 'border-green-400 bg-green-50' : 'border-gray-300 hover:border-indigo-400'}`}
            onDrop={(e) => handleDrop(e, setR2File)}
            onDragOver={handleDragOver}
            onClick={() => document.getElementById('r2-input')?.click()}
          >
            <input
              id="r2-input"
              type="file"
              className="hidden"
              accept=".fastq,.fq,.gz"
              onChange={(e) => e.target.files?.[0] && setR2File(e.target.files[0])}
            />
            {r2File ? (
              <>
                <Check className="mx-auto h-12 w-12 text-green-500" />
                <p className="mt-2 text-sm text-green-600">{r2File.name}</p>
              </>
            ) : (
              <>
                <File className="mx-auto h-12 w-12 text-gray-400" />
                <p className="mt-2 text-sm text-gray-600">R2 FASTQ (optional)</p>
                <p className="text-xs text-gray-500">Drag & drop or click to upload</p>
              </>
            )}
          </div>

          {/* Whitelist Upload */}
          <div
            className={`border-2 border-dashed rounded-lg p-6 text-center cursor-pointer
              ${whitelistFile ? 'border-green-400 bg-green-50' : 'border-gray-300 hover:border-indigo-400'}`}
            onDrop={(e) => handleDrop(e, setWhitelistFile)}
            onDragOver={handleDragOver}
            onClick={() => document.getElementById('whitelist-input')?.click()}
          >
            <input
              id="whitelist-input"
              type="file"
              className="hidden"
              accept=".txt,.tsv"
              onChange={(e) => e.target.files?.[0] && setWhitelistFile(e.target.files[0])}
            />
            {whitelistFile ? (
              <>
                <Check className="mx-auto h-12 w-12 text-green-500" />
                <p className="mt-2 text-sm text-green-600">{whitelistFile.name}</p>
              </>
            ) : (
              <>
                <File className="mx-auto h-12 w-12 text-gray-400" />
                <p className="mt-2 text-sm text-gray-600">Barcode Whitelist</p>
                <p className="text-xs text-gray-500">Drag & drop or click to upload</p>
              </>
            )}
          </div>
        </div>

        <div className="mt-6 flex justify-end">
          <button
            onClick={() => uploadMutation.mutate()}
            disabled={!r1File || uploadMutation.isPending}
            className="bg-indigo-600 text-white px-6 py-2 rounded-md hover:bg-indigo-700
              disabled:bg-gray-400 disabled:cursor-not-allowed flex items-center space-x-2"
          >
            {uploadMutation.isPending ? (
              <span>Uploading...</span>
            ) : (
              <>
                <UploadIcon className="h-4 w-4" />
                <span>Upload Files</span>
              </>
            )}
          </button>
        </div>

        {uploadMutation.isError && (
          <div className="mt-4 p-4 bg-red-50 rounded-md flex items-center space-x-2 text-red-700">
            <AlertCircle className="h-5 w-5" />
            <span>Upload failed. Please try again.</span>
          </div>
        )}

        {jobId && (
          <div className="mt-4 p-4 bg-green-50 rounded-md">
            <p className="text-green-700">
              Files uploaded successfully! Job ID: <code className="bg-green-100 px-1">{jobId}</code>
            </p>
            <p className="text-sm text-green-600 mt-2">
              Proceed to the Pipeline tab to start analysis.
            </p>
          </div>
        )}
      </div>
    </div>
  )
}
