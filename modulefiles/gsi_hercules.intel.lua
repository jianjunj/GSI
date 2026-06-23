help([[
]])

prepend_path("MODULEPATH", "/apps/contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.1.0/install/modulefiles/Core")

local stack_oneapi_ver=os.getenv("stack_oneapi_ver") or "2024.2.1"
local stack_intel_oneapi_mpi_ver=os.getenv("stack_intel_oneapi_mpi_ver") or "2021.13"
local intel_oneapi_mkl_ver=os.getenv("intel_oneapi_mkl_ver") or "2024.2.1"
local stack_python_ver=os.getenv("stack_python_ver") or "3.11.7"
local cmake_ver=os.getenv("cmake_ver") or "3.27.9"
local crtm_fix_ver=os.getenv("crtm_fix_ver") or "2.4.0.2"

load(pathJoin("stack-oneapi", stack_oneapi_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_intel_oneapi_mpi_ver))
load(pathJoin("intel-oneapi-mkl", intel_oneapi_mkl_ver))
load(pathJoin("stack-python", stack_python_ver))
load(pathJoin("cmake", cmake_ver))

load("gsi_common")
setenv("crtm_ROOT","/work/noaa/da/eliu/HERCULES/CRTM/crtm_v2.4.1-jedi.1/build")
setenv("crtm_VERSION","2.4.1-jedi.1")
setenv("CRTM_INC","/work/noaa/da/eliu/HERCULES/CRTM/crtm_v2.4.1-jedi.1/build/module")
setenv("CRTM_LIB","/work/noaa/da/eliu/HERCULES/CRTM/crtm_v2.4.1-jedi.1/build/lib/libcrtm_static.a")
setenv("CRTM_FIX","/work/noaa/da/eliu/HERCULES/CRTM-fix/crtm_v2.4.1-jedi.1-fix_gdasapp")
whatis("Name: crtm")
whatis("Version: 2.4.1-jedi.1")
whatis("Category: library")
whatis("Description: crtm library")

load(pathJoin("prod_util", prod_util_ver))
load("intel-oneapi-mkl/2022.2.1")

load("tar/1.34")

pushenv("CFLAGS", "-xHOST")
pushenv("FFLAGS", "-xHOST")

setenv("CC","mpiicc")
setenv("CXX","mpiicpc")
setenv("FC","mpiifort")

pushenv("GSI_BINARY_SOURCE_DIR", "/work2/noaa/global/role-global/fix/gsi/20251105")
setenv("CRTM_FIX", pathJoin("/work2/noaa/global/role-global/fix/crtm", "v" .. crtm_fix_ver))

whatis("Description: GSI environment on Hercules with Intel Compilers")
