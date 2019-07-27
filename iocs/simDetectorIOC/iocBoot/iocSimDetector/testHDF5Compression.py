# This program tests saving HDF5 files with compression and reading back the compressed files

import epics
import time
import h5py
import numpy as np
import os

CAM = '13SIM1:cam1:'
HDF5 = '13SIM1:HDF1:'
# For Linux
READ_PATH = '/home/epics/scratch/'
# For Windows
WRITE_PATH = 'J:/epics/scratch'

def write_read_file(filename, compression, bloscShuffle, bloscCompressor, bloscLevel, jpegQuality):
  t = epics.caput(HDF5 + 'FileName', filename)
  t = epics.caput(HDF5 + 'Compression', compression)
  print()
  print('Compression = ', compression)
  if (compression == 'Blosc'):
    t = epics.caput(HDF5 + 'BloscShuffle', bloscShuffle)
    t = epics.caput(HDF5 + 'BloscCompressor', bloscCompressor)
    t = epics.caput(HDF5 + 'BloscLevel', bloscLevel)
    print('BloscShuffle=', bloscShuffle, ' BloscCompressor=', bloscCompressor, ' BloscLevel=', bloscLevel)
  if (compression == 'JPEG'):
    t = epics.caput(HDF5 + 'JPEGQuality', jpegQuality)
    print('JPEG quality = ', jpegQuality)
  t = epics.caput(HDF5 + 'WriteFile', 1)
  time.sleep(0.5)
  file = READ_PATH + filename + '_001.h5'
  print('File size = ', os.path.getsize(file))
  hf = h5py.File(file, 'r')
  data = hf.get('/entry/data/data')
  data = np.array(data)
  print('Data[0-2] = ', data[0,0], data[1,0], data[2,0])
  hf.close()


t = epics.caput(CAM + "AcquireTime", 0.1)
t = epics.caput(CAM + "ImageMode", "Continuous")
t = epics.caput(CAM + "DataType", "UInt8")
t = epics.caput(CAM + "ColorMode", "Mono")
t = epics.caput(CAM + "Acquire", 1)

t = epics.caput(HDF5 + 'EnableCallbacks', 'Enable')
t = epics.caput(HDF5 + 'FilePath', WRITE_PATH)
t = epics.caput(HDF5 + 'FileNumber', 1)
t = epics.caput(HDF5 + 'AutoIncrement', 'No')
t = epics.caput(HDF5 + 'FileTemplate', '%s%s_%3.3d.h5')
t = epics.caput(HDF5 + 'FileWriteMode', 'Single')

write_read_file('Uncompressed', 'None', 0, 0, 0, 0)
write_read_file('Blosc_lz4_bit_5',   'Blosc', 'Bit',  'LZ4',  5, 0)
write_read_file('Blosc_zlib_byte_5', 'Blosc', 'Byte', 'ZLIB', 5, 0)
write_read_file('LZ4',               'LZ4',    0, 0, 0, 0)
write_read_file('BSLZ4',             'BSLZ4',  0, 0, 0, 0)
write_read_file('JPEG_Q90',          'JPEG',   0, 0, 0, 90)




