# Field API

Field API aims to ease the management and the transfer of data between CPUs and
GPUs for the Météo-France/ECMWF software.

The API is using fypp heavily to generate the code for several types and
dimensions. It might look complicated, but if you are just using the API then
you should not worry about it.

# Compilation

## Requirements
Building FIELD_API requires:
- A Fortran 2008 compliant compiler with support for:
  - OpenMP for CPU multi-threading
  - OpenACC for GPU offload (optional)
- CMake (>= 3.24)
- [ecbuild](https://github.com/ecmwf/ecbuild) (cloned if not found)
- [fypp](https://github.com/aradi/fypp) (cloned if not found)
- [fiat](https://github.com/ecmwf-ifs/fiat/) (optional)

To build FIELD_API without fiat, the path to the directory containing the utility modules `oml_mod.F90`, `abor1.F90` and `parkind1.F90` must be specified using the CMake variable `UTIL_MODULE_PATH`.

## Build and test
```
mkdir build
cd build
cmake .. # configure FIELD_API build
make # build FIELD_API
make install #optional, install FIELD_API
ctest #Optional, will run the tests
```

## Using FIELD_API
FIELD_API can be compiled as an external library and used in any CMake project simply by setting the environment variable `field_api_ROOT` to the FIELD_API builddir. Alternatively, FIELD_API can also be installed to a location on the `PATH`.

## Arch files
Architecture specific environment variables and compiler flags can be set by sourcing one of the `env.sh` files provided in the `arch` directory.

## Features
Features of FIELD_API can be toggled by passing the following argument to the CMake configure step: `-DENABLE_<FEATURE>=ON/OFF`. The table below lists all the available features:

| Feature | Default | Description |
|:--- |:--- |:--- |
| TESTS | ON | Build the testing suite. |
| ACC | ON | Enable the use of OpenACC for GPU offload. |

## Supported compilers
The library has been tested with the nvhpc toolkit from Nvidia, version 23.9
and is continually tested with newer releases. It has also been tested on CPU
(-DENABLE_ACC=OFF) with GCC 12 and Intel 2021. The CI testing (CPU-only for now) uses GNU 11.4.

# Field API types

Field API can encapsulate arrays of different types (REAL, INTEGER, LOGICAL)
and different dimensions (2D, 3D, 4D). Adding new dimensions would be trivial
and just a matter of editing the field\_definitions.fypp file. Adding the
encapsulation of new type would be easy if they are similar to the ones already
encapsulated.

Field API can be used to encapsulate data when used by a single thread or with
multiple threads. By default (persistent option unset or set to true), the last
dimension is used to store the thread number, and then each thread will access
only a single dimension of the array through a view.

```
SUBROUTINE SUB()
USE FIELD_MODULE
USE FIELD_FACTORY_MODULE
CLASS(FIELD_2RB), POINTER :: FO => NULL()
TYPE(FIELD_2D_VIEW_PTR) :: V

!Will create a field with the first dimension going from 1 to 10 and second from 1 to OMP_NUM_THREADS
CALL FIELD_NEW(FO, LBOUNDS=/1,1/, UBOUNDS=/10,1/)

DO IBLK=1,NBLKS
  V => FO%GET_VIEW(IBLK)
  !do stuff with v
ENDDO
  
CALL FIELD_DELETE(FO)
```

Furthermore field API provides two way of encapsulating the data: wrappers
and owners

## WRAPPER

The field API wrappers (eg. FIELD\_2D\_WRAPPER) provide a way to encapsulate
data which was already allocated before entering a part of the code which uses
field API. It is really just adding a wrapper around an array.

For instance, let say there are some data used in the code that has nothing to
do with GPUs. The data are allocated there and used there. But at some point,
maybe deeper in the code, those data are in fact needed on GPU. Instead of
having to declare a field API object high in the call stack, one could simply
declare a field api wrapper when needed.

```
SUBROUTINE SUB(MYDATA)
USE FIELD_MODULE
USE FIELD_FACTORY_MODULE
INTEGER, INTENT(INOUT) :: MYDATA(:,:)
CLASS(FIELD_2RB), POINTER :: FW => NULL()

!Wrap MYDATA into field wrapper FW
CALL FIELD_NEW(FW, DATA=MYDATA)

!do stuff

CALL FIELD_DELETE(FW)

!MYDATA is still accessible
MYDATA(1,2) = 7
```

## OWNER

The field API owners (eg. FIELD\_2D\_OWNER) provide a way to declare and
allocate  data to be used on CPU and GPU. The data is allocated by the API, the 
user doesn't have to allocate data by itself. Similarly, the user doesn't 
deallocate the data by himself, it is done by the API. When creating a owner,
the user will need to provide the two arrays used to specify the lower and
upper bounds of the array that will be created by field api.

```
SUBROUTINE SUB()
USE FIELD_MODULE
USE FIELD_FACTORY_MODULE
CLASS(FIELD_2RB), POINTER :: FO => NULL()

!Allocate data with field API 
!The allocated data will have a first dimension
!going from 1 to 10 and a second from 0 to 10.
CALL FIELD_NEW(FO, /1,0/, /10,10/, PERSISTENT=.FALSE.)

!do stuff

CALL FIELD_DELETE(FO)
!The data has now be freed on CPU and GPU and cannot be accessed anymore
```

### Delayed allocations

Field owners also provide a way to delay the allocation of data. The user can
ask a field owner to be created without allocating the data. The allocation
would then happen only if the data would be requested at some point, later in
the program. It can be useful if one doesn't want to waste memory on data that
might be only conditionally used. But please keep in mind, that allocating data
can be slow and will slow down the program if done during a computation heavy
part of the code.

```
SUBROUTINE SUB(MYTEST)
LOGICAL, INTENT(IN) :: MYTEST
USE FIELD_MODULE
USE FIELD_FACTORY_MODULE
CLASS(FIELD_2RB), POINTER :: FO => NULL()

!Declare a field owner, no allocation will happen here
CALL FIELD_NEW(FO, /1,0/, /10,10/, PERSISTENT=.FALSE., DELAYED=.TRUE.)

IF (MYTEST) THEN
!do stuff with FO
!allocation wil happen here
ENDIF

CALL FIELD_DELETE(FO)
!The data will be freed if MYTEST was true, otherwise there are no data to deallocate
```

### Initialisation

In the case of field owner it is possible to initiliase it with a specific
value at creation time by adding the INIT\_VALUE optional argument.

```
   CLASS(FIELD_2IM), POINTER :: O => NULL()
   !This field owner will be initialised to 3
   CALL FIELD_NEW(O, LBOUNDS=[1,1], UBOUNDS=[10,10], INIT_VALUE=3_JPIM)
```

It is also possible to activate a debug value to initialise all non-initialised
owner. To do so it is necessary to import the module *field_init_debug_module*
and set *use_init_debug_value* to true. Then one can *set
init_debug_value_jpim* to a custom value.

```
   USE FIELD_INIT_DEBUG_VALUE_MODULE
   USE_INIT_DEBUG_VALUE = .TRUE.
   INIT_DEBUG_VALUE_JPIM = -7
   !This field owner will be initialised to -7
   CALL FIELD_NEW(O, LBOUNDS=[1,1], UBOUNDS=[10,10])
```

## Asynchronism

This functionnality is still being tested.

By default all data transfers are synchronous. So every call to the subroutines
GET\_HOST\_DATA, GET\_DEVICE\_DATA, SYNC\_HOST, SYNC\_DEVICE will stop the
program until the data are actually transfered. But sometimes it is possible to
interleave the data transfer with the computations. To do so you can add the
QUEUE parameter when calling the aforementioned subroutines. With this QUEUE
parameter the user will specify on which queue he wants the data transfer to
happen, and the subroutines will return without waiting for the data transfer
to finish. It is up to the user to be sure the data transfer has been done when
he actually wants to use the data. This can be checked by using the
WAIT\_FOR\_ASYNC\_QUEUE subroutine.

```
SUBROUTINE SUB(MYTEST)
USE FIELD_MODULE
USE FIELD_FACTORY_MODULE
CLASS(FIELD_2RB), POINTER :: FO => NULL()
CLASS(FIELD_2RB), POINTER :: FO2 => NULL()
LOGICAL, INTENT(IN) :: MYTEST

CALL FIELD_NEW(FO, /1,0/, /10,10/)
CALL FIELD_NEW(FO2, /1,0/, /10,10/)

!Do stuff with FO on GPUs
!Then transfer data to CPU
CALL FO%SYNC_HOST_RDONLY(QUEUE=2)

!Do stuff with FO2 on GPUs
!We didn't have to wait for the data transfer of FO to finish

!Make sure the data transfer for FO is finished
CALL WAIT_FOR_ASYNC_QUEUE(QUEUE=2)

...

```

## Statistics

Each field API variable maintains statistics about the time it spend on data
transfer and the number of time it happened. You can access them through the
field ```FW%STATS```

For instance:
```
...
NUM_CPU_GPU_TR=FW%STATS%TRANSFER_CPU_TO_GPU
AVG=FW%STATS%TOTAL_TIME_TRANSFER_CPU_TO_GPU/FW%STATS%TRANSFER_CPU_TO_GPU
write(*,*)"Num transfer CPU->GPU", NUM_CPU_GPU_TR
write(*,*)"Total/Avg Time spend on transfer CPU->GPU", NUM_CPU_GPU_TR, "/" AVG, 
...
```

## Note on GET\_VIEW

GET\_VIEW must only be called in sections of code running on the host. The
field's data must be present on the host. It will not work if the data are on
the device or if the field has not been allocated yet (when using the DELAY
option).

# Public API

For field api type:
```
SUBROUTINE FIELD_NEW(SELF, ...)
SUBROUTINE FIELD_DELETE(SELF)
SUBROUTINE DELETE_DEVICE
FUNCTION GET_VIEW(SELF, BLOCK_INDEX, ZERO) RESULT(VIEW_PTR)
SUBROUTINE GET_DEVICE_DATA_RDONLY (SELF, PPTR, QUEUE)
SUBROUTINE GET_DEVICE_DATA_RDWR (SELF, PPTR, QUEUE)
SUBROUTINE GET_HOST_DATA_RDONLY (SELF, PPTR, QUEUE)
SUBROUTINE GET_HOST_DATA_RDWR (SELF, PPTR, QUEUE)
SUBROUTINE SYNC_HOST_RDWR (SELF, QUEUE)
SUBROUTINE SYNC_HOST_RDONLY (SELF, QUEUE)
SUBROUTINE SYNC_DEVICE_RDWR (SELF, QUEUE)
SUBROUTINE SYNC_DEVICE_RDONLY (SELF, QUEUE)
SUBROUTINE COPY_OBJECT (SELF, LDCREATED)
SUBROUTINE WIPE_OBJECT (SELF, LDDELETED)
```

Utils:

```
SUBROUTINE WAIT_FOR_ASYNC_QUEUE(QUEUE)
TYPE FIELD_*D_PTR
```

Stats:
```
INTEGER :: TRANSFER_CPU_TO_GPU 
INTEGER :: TRANSFER_GPU_TO_CPU
REAL :: TOTAL_TIME_TRANSFER_CPU_TO_GPU
REAL :: TOTAL_TIME_TRANSFER_GPU_TO_CPU
```



