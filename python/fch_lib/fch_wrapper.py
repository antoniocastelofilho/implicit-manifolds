# fch_wrapper.py
import os
import ctypes as C
import numpy as np
import platform
from typing import Callable, Optional
import psutil

from numba import cfunc, njit, types, carray


class FCH_CPU_Runner:
    """
    Python wrapper for the C++ FCH CPU library.

    It loads the shared library with hardcoded paths, preloads dependencies,
    sets ctypes signatures, and compiles a GIL-free C callback via Numba so 
    your C++ threadpool can call it in parallel.

    You can pass either:
      - a plain Python math function: def f(x: np.ndarray, y: np.ndarray) -> None
        (must be Numba-nopython friendly; use math.*, not numpy.* inside), or
      - a precompiled Numba cfunc object (has `.address`).
    """

    # ---- C types ----
    _DoubleP = C.POINTER(C.c_double)
    _UInt32P = C.POINTER(C.c_uint32)
    _CALLBACK_VOIDP = C.c_void_p

    # ---- Cross-platform library paths ----
    @staticmethod
    def _get_library_path(build_path):
        """Get the appropriate library path based on the current platform."""
        system = platform.system().lower()
        
        if system == "linux":
            return build_path + "/liblib_fch_cpu.so"
        elif system == "windows":
            return build_path + "/lib_fch_cpu.dll"  # Note: no 'lib' prefix on Windows
        elif system == "darwin":  # macOS
            return build_path + "/liblib_fch_cpu.dylib"
        else:
            raise OSError(f"Unsupported platform: {system}")

    def __init__(self, build_path = "./build/"):
        # Load main library with platform-appropriate extension
        library_path = self._get_library_path(build_path)
        
        if not os.path.exists(library_path):
            raise FileNotFoundError(f"Main library not found at: {library_path}")
        
        self._lib = C.CDLL(library_path)

        # Get optional symbols safely
        # FCH with threadpool
        self._c_threadpool_func = getattr(self._lib, "fch_cpu_threadpooljobs", None)
        if self._c_threadpool_func is not None:
            self._setup_threadpool_signature()

        # FCH batched
        self._c_batched_func = getattr(self._lib, "fch_cpu_batched", None)
        if self._c_batched_func is not None:
            self._setup_batched_signature()

        # Combinatorial skeleton
        self._c_skeleton = getattr(self._lib, "skeleton", None)
        if self._c_skeleton is not None:
            self._setup_skeleton_signature()

        # Simplicial decomposition
        self._c_skeleton_simplicial_decomp = getattr(self._lib, "skeleton_simplicial_decomp", None)
        if self._c_skeleton_simplicial_decomp is not None:
            self._setup_skeleton_simplicial_decomp_signature()

        # Find first simplex
        self._c_find_first_simplex = getattr(self._lib, "find_first_simplex", None)
        if self._c_find_first_simplex is not None:
            self._setup_find_first_simplex_signature()

        # Get coface memory byte size
        self._c_get_coface_block_alloc_size = getattr(self._lib, "get_coface_block_alloc_size", None)
        if self._c_get_coface_block_alloc_size is not None:
            self._setup_get_coface_block_alloc_size_signature()

        if (self._c_threadpool_func is None) or (self._c_batched_func is None) or (self._c_skeleton is None) or (self._c_skeleton_simplicial_decomp is None) or (self._c_find_first_simplex is None) or (self._c_get_coface_block_alloc_size is None):
            raise AttributeError(
                "Library does not export 'fch_cpu_threadpooljobs', 'fch_cpu_batched', 'skeleton', 'skeleton_simplicial_decomp', 'find_first_simplex', or 'get_coface_block_alloc_size'."
            )

        # Keep compiled callbacks alive across calls (optional cache)
        self._keepalives = []


    # ----------------- ctypes signatures -----------------

    def _setup_threadpool_signature(self):
        if self._c_threadpool_func is None:
            return
        self._c_threadpool_func.argtypes = [
            C.c_uint32, C.c_uint32,                # NDIM, KDIM
            self._DoubleP, self._DoubleP, self._DoubleP, self._DoubleP,  # DOMAIN_{MIN,MIN_CLIP,MAX,MAX_CLIP}
            self._UInt32P,                         # DOMAIN_DIV
            self._CALLBACK_VOIDP,                  # func (raw function pointer)
            self._UInt32P,                         # FIRST_COORD
            self._UInt32P,                         # FIRST_EIDX
            C.c_uint32,                            # n_threads
            C.c_char_p,                            # filepath_out
            C.c_uint32,                            # print_interval
        ]
        self._c_threadpool_func.restype = None

    def _setup_batched_signature(self):
        if self._c_batched_func is None:
            return
        self._c_batched_func.argtypes = [
            C.c_uint32, C.c_uint32,                # NDIM, KDIM
            self._DoubleP, self._DoubleP, self._DoubleP, self._DoubleP,  # DOMAIN_{MIN,MIN_CLIP,MAX,MAX_CLIP}
            self._UInt32P,                         # DOMAIN_DIV
            self._CALLBACK_VOIDP,                  # func (raw function pointer)
            self._UInt32P,                         # FIRST_COORD
            self._UInt32P,                         # FIRST_EIDX
            C.c_uint32,                            # n_threads
            C.c_uint32,                            # n_cofaces
            C.c_char_p,                            # filepath_out
        ]
        self._c_batched_func.restype = None

    # NEW: skeleton signature
    def _setup_skeleton_signature(self):
        if self._c_skeleton is None:
            return
        self._c_skeleton.argtypes = [
            C.c_char_p,  # filepath_in
            C.c_char_p,  # filepath_out
        ]
        self._c_skeleton.restype = None  # adjust if your C++ returns something else

    def _setup_skeleton_simplicial_decomp_signature(self):
        if self._c_skeleton_simplicial_decomp is None:
            return
        self._c_skeleton_simplicial_decomp.argtypes = [
            C.c_char_p,  # filepath_in
            C.c_char_p,  # filepath_out
            C.c_uint8,   # is_output_binary
        ]
        self._c_skeleton_simplicial_decomp.restype = None  # adjust if your C++ returns something else

    def _setup_find_first_simplex_signature(self):
        if self._c_find_first_simplex is None:
            return
        self._c_find_first_simplex.argtypes = [
            C.c_uint32, C.c_uint32,                # NDIM, KDIM
            self._DoubleP, self._DoubleP,          # DOMAIN_{MIN,MAX}
            self._UInt32P,                         # DOMAIN_DIV
            self._CALLBACK_VOIDP,                  # func (raw function pointer)
            self._DoubleP,                         # PRE_ROOT
            C.c_uint32,                            # PRE_ROOT_LENGTH
            C.c_uint32,                            # n_threads
            C.c_char_p,                            # filepath_out
            C.c_uint32,                            # initial_exponent
            C.c_uint32,                            # max_exponent
            C.c_uint32,                            # max_iter_per_exponent
        ]
        self._c_find_first_simplex.restype = C.c_int32

    def _setup_get_coface_block_alloc_size_signature(self):
        if self._c_get_coface_block_alloc_size is None:
            return
        self._c_get_coface_block_alloc_size.argtypes = [
            C.c_uint32, C.c_uint32,                # NDIM, KDIM
        ]
        self._c_get_coface_block_alloc_size.restype = C.c_uint64

    # ----------------- callback compilation -----------------

    @staticmethod
    def _as_c_contig(a: np.ndarray, dtype) -> np.ndarray:
        return np.ascontiguousarray(a, dtype=dtype)

    @staticmethod
    def _compile_user_math(user_math_func: Callable):
        """
        Compile user's math (x: float64[:], y: float64[:]) to nopython.
        User code must use math.* (not numpy.*) inside for best compatibility.
        """
        sig = types.void(types.float64[:], types.float64[:])
        return njit(sig, nogil=True, fastmath=True, cache=True)(user_math_func)

    @staticmethod
    def _build_numba_cfunc_from_njit(njit_func, NDIM: int, KDIM: int):
        """
        Wrap an njit function into a true C callback (no GIL).
        """
        sig_c = types.void(types.CPointer(types.double), types.CPointer(types.double))

        @cfunc(sig_c, nopython=True, fastmath=True, cache=True)
        def _shim(x_ptr, y_ptr):
            x = carray(x_ptr, NDIM)
            y = carray(y_ptr, KDIM)
            njit_func(x, y)

        cfunc_obj = _shim # to prevent pylance from complaining: assigning _shim to an object, and then calling the address attribute of the object.
        return _shim, C.c_void_p(cfunc_obj.address)

    @staticmethod
    def _accept_user_callback(math_func_or_cfunc, NDIM: int, KDIM: int):
        """
        Accept either:
          - a Numba cfunc object (has .address), or
          - a plain Python function (compile to njit, then wrap in cfunc).
        Returns (keepalive_obj, c_void_p_address).
        """
        # Case 1: already a Numba cfunc
        if hasattr(math_func_or_cfunc, "address"):
            cfun = math_func_or_cfunc
            return cfun, C.c_void_p(cfun.address)

        # Case 2: compile Python function
        njitted = FCH_CPU_Runner._compile_user_math(math_func_or_cfunc)
        return FCH_CPU_Runner._build_numba_cfunc_from_njit(njitted, NDIM, KDIM)

    # ----------------- public API -----------------

    def run_threadpool(
        self,
        *,
        math_func: Callable,                   # See _accept_user_callback docs
        domain_min: np.ndarray,
        domain_max: np.ndarray,
        domain_div: np.ndarray,
        first_coord: np.ndarray,
        first_eidx: np.ndarray,
        out_path: str,
        n_threads: int = 0,
        domain_min_clip: Optional[np.ndarray] = None,
        domain_max_clip: Optional[np.ndarray] = None,
        print_interval: int = 0
    ):
        """Call fch_cpu_threadpooljobs with a GIL-free callback."""
        if self._c_threadpool_func is None:
            raise NotImplementedError("'fch_cpu_threadpooljobs' not in library.")

        NDIM = int(domain_min.shape[0])
        KDIM = int(first_eidx.shape[0])

        if n_threads < 1:
            n_threads = psutil.cpu_count(logical=False) or 1
            print(f"n_threads set automatically: {n_threads}")

        if domain_min_clip is None:
            domain_min_clip = domain_min
        if domain_max_clip is None:
            domain_max_clip = domain_max

        # Shape checks
        assert domain_max.shape[0]      == NDIM, "domain_max shape mismatch"
        assert domain_min_clip.shape[0] == NDIM, "domain_min_clip shape mismatch"
        assert domain_max_clip.shape[0] == NDIM, "domain_max_clip shape mismatch"
        assert domain_div.shape[0]      == NDIM, "domain_div shape mismatch"
        assert first_coord.shape[0]     == NDIM, "first_coord shape mismatch"

        # Contiguous inputs
        dmin      = self._as_c_contig(domain_min,      np.float64)
        dmin_clip = self._as_c_contig(domain_min_clip, np.float64)
        dmax      = self._as_c_contig(domain_max,      np.float64)
        dmax_clip = self._as_c_contig(domain_max_clip, np.float64)
        ddiv      = self._as_c_contig(domain_div,      np.uint32)
        fcoord    = self._as_c_contig(first_coord,     np.uint32)
        feidx     = self._as_c_contig(first_eidx,      np.uint32)
        out_bytes = os.fspath(out_path).encode("utf-8")

        # Build (or accept) a GIL-free callback
        keepalive, cb_addr = self._accept_user_callback(math_func, NDIM, KDIM)
        self._keepalives.append(keepalive)

        # Call
        self._c_threadpool_func(
            C.c_uint32(NDIM), C.c_uint32(KDIM),
            dmin.ctypes.data_as(self._DoubleP),
            dmin_clip.ctypes.data_as(self._DoubleP),
            dmax.ctypes.data_as(self._DoubleP),
            dmax_clip.ctypes.data_as(self._DoubleP),
            ddiv.ctypes.data_as(self._UInt32P),
            cb_addr,
            fcoord.ctypes.data_as(self._UInt32P),
            feidx.ctypes.data_as(self._UInt32P),
            C.c_uint32(n_threads),
            out_bytes,
            C.c_uint32(print_interval),
        )

    def run_batched(
        self,
        *,
        math_func: Callable,                   # See _accept_user_callback docs
        domain_min: np.ndarray,
        domain_max: np.ndarray,
        domain_div: np.ndarray,
        first_coord: np.ndarray,
        first_eidx: np.ndarray,
        out_path: str,
        n_cofaces: int = 0,
        n_threads: int = 0,
        domain_min_clip: Optional[np.ndarray] = None,
        domain_max_clip: Optional[np.ndarray] = None,
    ):
        """Call fch_cpu_batched with a GIL-free callback."""
        if self._c_batched_func is None:
            raise NotImplementedError("'fch_cpu_batched' not in library.")

        NDIM = int(domain_min.shape[0])
        KDIM = int(first_eidx.shape[0])

        if n_threads < 1:
            n_threads = psutil.cpu_count(logical=False) or 1
            print(f"n_threads set automatically: {n_threads}")

        if n_cofaces < 1:
            # compute the total amount of cofaces that could be created with the available memory.
            available_memory_bytes = psutil.virtual_memory().available
            coface_block_alloc_size = self.get_coface_block_alloc_size(NDIM, KDIM)
            # print(available_memory_bytes)
            # print(coface_block_alloc_size)
            n_cofaces = available_memory_bytes // coface_block_alloc_size
            if n_cofaces < 2:
                raise ValueError("not enough memory")
            # use only half of the available memory for the cofaces, letting the other half to be used by the hash table.
            n_cofaces = max(1, n_cofaces // 2)
            print(f"n_cofaces set automatically: {n_cofaces}")


        if domain_min_clip is None:
            domain_min_clip = domain_min
        if domain_max_clip is None:
            domain_max_clip = domain_max

        # Shape checks
        assert domain_max.shape[0]      == NDIM, "domain_max shape mismatch"
        assert domain_min_clip.shape[0] == NDIM, "domain_min_clip shape mismatch"
        assert domain_max_clip.shape[0] == NDIM, "domain_max_clip shape mismatch"
        assert domain_div.shape[0]      == NDIM, "domain_div shape mismatch"
        assert first_coord.shape[0]     == NDIM, "first_coord shape mismatch"

        # Contiguous inputs
        dmin      = self._as_c_contig(domain_min,      np.float64)
        dmin_clip = self._as_c_contig(domain_min_clip, np.float64)
        dmax      = self._as_c_contig(domain_max,      np.float64)
        dmax_clip = self._as_c_contig(domain_max_clip, np.float64)
        ddiv      = self._as_c_contig(domain_div,      np.uint32)
        fcoord    = self._as_c_contig(first_coord,     np.uint32)
        feidx     = self._as_c_contig(first_eidx,      np.uint32)
        out_bytes = os.fspath(out_path).encode("utf-8")

        # Build (or accept) a GIL-free callback
        keepalive, cb_addr = self._accept_user_callback(math_func, NDIM, KDIM)
        self._keepalives.append(keepalive)

        # Call
        self._c_batched_func(
            C.c_uint32(NDIM), C.c_uint32(KDIM),
            dmin.ctypes.data_as(self._DoubleP),
            dmin_clip.ctypes.data_as(self._DoubleP),
            dmax.ctypes.data_as(self._DoubleP),
            dmax_clip.ctypes.data_as(self._DoubleP),
            ddiv.ctypes.data_as(self._UInt32P),
            cb_addr,
            fcoord.ctypes.data_as(self._UInt32P),
            feidx.ctypes.data_as(self._UInt32P),
            C.c_uint32(n_threads),
            C.c_uint32(n_cofaces),
            out_bytes,
        )

    # ----------------- NEW public API: skeleton -----------------

    def skeleton(self, *, in_path: str, out_path: str):
        """
        Call `skeleton(const char* filepath_in, const char* filepath_out)` from the C++ library.

        Parameters
        ----------
        in_path : str
            Input file path expected by the C++ implementation.
        out_path : str
            Output file path to be produced by the C++ implementation.
        """
        if self._c_skeleton is None:
            raise NotImplementedError("'skeleton' not found in the loaded library.")

        # Encode as UTF-8 bytes for c_char_p
        in_bytes = os.fspath(in_path).encode("utf-8")
        out_bytes = os.fspath(out_path).encode("utf-8")

        self._c_skeleton(in_bytes, out_bytes)

    def skeleton_simplicial_decomp(self, *, in_path: str, out_path: str, is_output_binary: bool = False):
        """
        Call `skeleton_simplicial_decomp(const char* filepath_in, const char* filepath_out)` from the C++ library.

        Parameters
        ----------
        in_path : str
            Input file path expected by the C++ implementation.
        out_path : str
            Output file path to be produced by the C++ implementation.
        """
        if self._c_skeleton_simplicial_decomp is None:
            raise NotImplementedError("'skeleton_simplicial_decomp' not found in the loaded library.")

        # Encode as UTF-8 bytes for c_char_p
        in_bytes = os.fspath(in_path).encode("utf-8")
        out_bytes = os.fspath(out_path).encode("utf-8")

        is_output_binary_uint8 = C.c_uint8(0)
        if is_output_binary:
            is_output_binary_uint8 = C.c_uint8(1)

        self._c_skeleton_simplicial_decomp(in_bytes, out_bytes, is_output_binary_uint8)

    def find_first_simplex(
        self,
        *,
        ndim: int,
        kdim: int,
        math_func: Callable,                   # See _accept_user_callback docs
        domain_min: np.ndarray,
        domain_max: np.ndarray,
        domain_div: np.ndarray,
        out_path: str,
        n_threads: int = 0,
        PRE_ROOT: Optional[np.ndarray] = None,
        PRE_ROOT_LENGTH: int = -1,
        initial_exponent: int = 4,
        max_exponent: int = 8,
        max_iter_per_exponent: int = 20000
    ):
        """Call fch_cpu_threadpooljobs with a GIL-free callback."""
        if self._c_find_first_simplex is None:
            raise NotImplementedError("'find_first_simplex' not in library.")

        NDIM = int(domain_min.shape[0])
        if ndim != NDIM:
            raise ValueError("ndim != NDIM")
        
        KDIM = kdim

        if n_threads < 1:
            n_threads = psutil.cpu_count(logical=False) or 1
            print(f"n_threads set automatically: {n_threads}")

        # Shape checks
        assert domain_min.shape[0]      == NDIM, "domain_min shape mismatch"
        assert domain_max.shape[0]      == NDIM, "domain_max shape mismatch"
        assert domain_div.shape[0]      == NDIM, "domain_div shape mismatch"
        if (PRE_ROOT is not None) and (PRE_ROOT_LENGTH >= 0):
            assert PRE_ROOT.shape[0]    >= PRE_ROOT_LENGTH, "PRE_ROOT shape mismatch"

        # Contiguous inputs
        dmin      = self._as_c_contig(domain_min,      np.float64)
        dmax      = self._as_c_contig(domain_max,      np.float64)
        ddiv      = self._as_c_contig(domain_div,      np.uint32)
        root0 = None
        proot0 = None
        if PRE_ROOT is not None:
            root0 = self._as_c_contig(PRE_ROOT,      np.float64)
            proot0 = root0.ctypes.data_as(self._DoubleP)
            if PRE_ROOT_LENGTH < 0:
                PRE_ROOT_LENGTH = PRE_ROOT.shape[0]
        else:
            PRE_ROOT_LENGTH = 0

        out_bytes = os.fspath(out_path).encode("utf-8")

        # Build (or accept) a GIL-free callback
        keepalive, cb_addr = self._accept_user_callback(math_func, NDIM, KDIM)
        self._keepalives.append(keepalive)

        # Call
        return self._c_find_first_simplex(
            C.c_uint32(NDIM), 
            C.c_uint32(KDIM),
            dmin.ctypes.data_as(self._DoubleP),
            dmax.ctypes.data_as(self._DoubleP),
            ddiv.ctypes.data_as(self._UInt32P),
            cb_addr,
            proot0,
            C.c_uint32(PRE_ROOT_LENGTH),
            C.c_uint32(n_threads),
            out_bytes,
            C.c_uint32(initial_exponent),
            C.c_uint32(max_exponent),
            C.c_uint32(max_iter_per_exponent)
        )
    
    def read_find_first_simplex_binary(self,  filepath: str):
        """
        Read the binary file written by the C++ find_first_simplex() implementation.

        Layout (little-endian):
        uint32 ndim
        uint32 kdim
        double[ndim]  dmin_shifted   = dmin_fin + offset
        double[ndim]  dmin_orig      = original dmin (clip)
        double[ndim]  dmax_shifted   = dmax_fin + offset
        double[ndim]  dmax_orig      = original dmax (clip)
        uint32[ndim]  ddiv
        uint32[ndim]  grid_coord
        uint32[kdim]  eidx

        Returns
        -------
        dict with:
        - 'ndim' : int
        - 'kdim' : int
        - 'dmin_shifted' : np.ndarray, dtype=float64, shape (ndim,)
        - 'dmin_orig'    : np.ndarray, dtype=float64, shape (ndim,)
        - 'dmax_shifted' : np.ndarray, dtype=float64, shape (ndim,)
        - 'dmax_orig'    : np.ndarray, dtype=float64, shape (ndim,)
        - 'ddiv'         : np.ndarray, dtype=uint32,  shape (ndim,)
        - 'grid_coord'   : np.ndarray, dtype=uint32,  shape (ndim,)
        - 'eidx'         : np.ndarray, dtype=uint32,  shape (kdim,)
        """
        with open(filepath, "rb") as f:
            buf = f.read()

        mv = memoryview(buf)
        pos = 0

        def read_u32():
            nonlocal pos
            if pos + 4 > len(mv):
                raise ValueError("Unexpected EOF while reading uint32")
            val = int.from_bytes(mv[pos:pos+4], "little", signed=False)
            pos += 4
            return val

        def read_f64_array(n):
            nonlocal pos
            nbytes = n * 8
            if pos + nbytes > len(mv):
                raise ValueError("Unexpected EOF while reading float64 array")
            arr = np.frombuffer(mv[pos:pos+nbytes], dtype="<f8", count=n)
            pos += nbytes
            # copy so the result isn't tied to the underlying memoryview
            return arr.copy()

        def read_u32_array(n):
            nonlocal pos
            nbytes = n * 4
            if pos + nbytes > len(mv):
                raise ValueError("Unexpected EOF while reading uint32 array")
            arr = np.frombuffer(mv[pos:pos+nbytes], dtype="<u4", count=n)
            pos += nbytes
            return arr.copy()

        ndim = read_u32()
        kdim = read_u32()

        dmin_shifted = read_f64_array(ndim)
        dmin_orig    = read_f64_array(ndim)
        dmax_shifted = read_f64_array(ndim)
        dmax_orig    = read_f64_array(ndim)
        ddiv         = read_u32_array(ndim)
        grid_coord   = read_u32_array(ndim)
        eidx         = read_u32_array(kdim)

        # Optional strict size check:
        expected_size = (2*4) + (4*ndim*8) + (2*ndim*4) + (kdim*4)
        if len(mv) != expected_size:
            # Not fatal, but you can make this a hard error if desired
            # raise ValueError(f"File size mismatch: got {len(mv)} bytes, expected {expected_size}")
            pass

        return {
            "ndim": int(ndim),
            "kdim": int(kdim),
            "dmin":         dmin_shifted,
            "dmin_clip":    dmin_orig,
            "dmax":         dmax_shifted,
            "dmax_clip":    dmax_orig,
            "ddiv":         ddiv,
            "first_coord":  grid_coord,
            "first_eidx":   eidx,
        }
    
    def read_fch_output_binary(self, filepath: str) -> np.ndarray:
        """
        Read a binary file laid out as:
        Header (little-endian uint32 fields):
            ndim, kdim, total_count, offset0, offset1, offset2, offset3, offset4, block_size
        Then:
            total_count blocks, each block_size bytes.

        Within each block, from byte offset "offset3" (relative to the block start)
        to the end of the block there are exactly (2 * ndim) doubles (float64).
        This function extracts those doubles from each block in a strided fashion
        and returns them as a (total_count, 2*ndim) array of dtype float64.

        Parameters
        ----------
        filepath : str
            Path to the binary file.

        Returns
        -------
        np.ndarray
            Array of shape (total_count, 2*ndim), dtype float64.
        """
        with open(filepath, "rb") as f:
            buf = f.read()

        mv = memoryview(buf)
        pos = 0

        def read_u32():
            nonlocal pos
            if pos + 4 > len(mv):
                raise ValueError("Unexpected EOF while reading uint32")
            val = int.from_bytes(mv[pos:pos+4], "little", signed=False)
            pos += 4
            return val

        # ---- parse header ----
        ndim = read_u32()
        kdim = read_u32()             # parsed but not used here
        total_count = read_u32()
        div = []
        for i in range(ndim):
            div.append(read_u32())
        offset0 = read_u32()
        offset1 = read_u32()
        offset2 = read_u32()
        offset3 = read_u32()
        offset4 = read_u32()
        block_size = read_u32()

        header_size = (9 + ndim) * 4  # 9 + ndim uint32 fields

        # ---- basic sanity checks ----
        if offset3 >= block_size:
            raise ValueError(f"offset3 ({offset3}) must be < block_size ({block_size})")

        row_bytes = (2 * ndim) * 8  # bytes for 2*ndim float64s
        if offset3 + row_bytes > block_size:
            raise ValueError(
                f"Block payload too small: offset3 ({offset3}) + {row_bytes} bytes "
                f"exceeds block_size ({block_size})."
            )

        total_needed = header_size + total_count * block_size
        if len(mv) < total_needed:
            raise ValueError(
                f"File too small: {len(mv)} bytes; need at least {total_needed} bytes "
                f"for {total_count} blocks of size {block_size}."
            )

        # ---- read rows in a strided manner ----
        out = np.empty((total_count, 2 * ndim), dtype=np.float64)
        for i in range(total_count):
            block_start = header_size + i * block_size
            data_start = block_start + offset3
            # Create a view on just this row's bytes, then decode as little-endian float64
            row_mv = mv[data_start : data_start + row_bytes]
            out[i] = np.frombuffer(row_mv, dtype="<f8", count=2 * ndim)

        return out
    
    def get_coface_block_alloc_size(
        self,
        ndim: int,
        kdim: int):
        
        if self._c_get_coface_block_alloc_size is None:
            raise NotImplementedError("'get_coface_block_alloc_size' not found in the loaded library.")

        return self._c_get_coface_block_alloc_size(
            C.c_uint32(ndim), 
            C.c_uint32(kdim),
            )

    def read_simplicial_decomp_binary(self, filepath: str):
        """
        Parse a binary file with layout:
        uint32  ndim
        uint32  kdim
        uint32[ndim]                  domain_div
        uint32  hyper_count
        uint32[hyper_count * ndim]    hypercubes (row-major: hyper_count rows, ndim cols)
        uint32  vertex_count
        double[vertex_count * ndim]   vertices  (row-major: vertex_count rows, ndim cols)
        uint32  simplex_count
        int64[simplex_count * L]      simplices (row-major: simplex_count rows, L cols)
            where L = 2*(ndim - kdim + 1) + 2

        Returns
        -------
        domain_div : np.ndarray, dtype=np.uint32, shape (ndim,)
        hypercubes : np.ndarray, dtype=np.uint32, shape (hyper_count, ndim)
        vertices   : np.ndarray, dtype=np.float64, shape (vertex_count, ndim)
        simplices  : np.ndarray, dtype=np.int64,   shape (simplex_count, L)
        """
        with open(filepath, "rb") as f:
            buf = f.read()

        mv = memoryview(buf)
        pos = 0

        def need(nbytes: int):
            if pos + nbytes > len(mv):
                raise ValueError(f"Unexpected EOF: need {nbytes} more bytes at position {pos}")

        def read_u32() -> int:
            nonlocal pos
            need(4)
            val = int.from_bytes(mv[pos:pos+4], "little", signed=False)
            pos += 4
            return val

        def read_u32_array(n: int) -> np.ndarray:
            nonlocal pos
            nbytes = 4 * n
            need(nbytes)
            arr = np.frombuffer(mv[pos:pos+nbytes], dtype="<u4", count=n)
            pos += nbytes
            return arr.copy()

        def read_f64_array(n: int) -> np.ndarray:
            nonlocal pos
            nbytes = 8 * n
            need(nbytes)
            arr = np.frombuffer(mv[pos:pos+nbytes], dtype="<f8", count=n)
            pos += nbytes
            return arr.copy()

        def read_i64_array(n: int) -> np.ndarray:
            nonlocal pos
            nbytes = 8 * n
            need(nbytes)
            arr = np.frombuffer(mv[pos:pos+nbytes], dtype="<i8", count=n)
            pos += nbytes
            return arr.copy()

        # ---- header ----
        ndim = read_u32()
        kdim = read_u32()

        # ---- domain_div ----
        domain_div = read_u32_array(ndim)

        # ---- hypercubes ----
        hyper_count = read_u32()
        if hyper_count < 0:
            raise ValueError("hyper_count must be non-negative")
        hc_flat = read_u32_array(hyper_count * ndim)
        hypercubes = hc_flat.reshape(hyper_count, ndim)

        # ---- vertices ----
        vertex_count = read_u32()
        if vertex_count < 0:
            raise ValueError("vertex_count must be non-negative")
        vtx_flat = read_f64_array(vertex_count * ndim)
        vertices = vtx_flat.reshape(vertex_count, ndim)

        # ---- simplices ----
        simplex_count = read_u32()
        if simplex_count < 0:
            raise ValueError("simplex_count must be non-negative")
        L = 2 * (ndim - kdim + 1) + 2
        if L <= 0:
            raise ValueError(f"Computed simplex row length L={L} must be positive (ndim={ndim}, kdim={kdim})")
        simp_flat = read_i64_array(simplex_count * L)
        simplices = simp_flat.reshape(simplex_count, L)

        # Optional: ensure we consumed the whole file (comment out if not desired)
        # if pos != len(mv):
        #     raise ValueError(f"Trailing data: consumed {pos} bytes, file has {len(mv)} bytes")

        return ndim, kdim, domain_div, hypercubes, vertices, simplices