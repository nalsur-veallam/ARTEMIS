import ctypes

test = ctypes.CDLL('denoise')

test.main()
