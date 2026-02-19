/*
 * (C) Copyright 2014- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/*!
 * @mainpage
 * This project declares the C-API to the IFS trans-library.\n
 * For documentation of all available functions, see @ref trans.h .
 *
 * @section About
 *
 * This library gives access to spectral transforms on the sphere.
 * The library is capable to take advantage of a MPI-distributed-memory environment,
 * and can use OpenMP-shared-memory parallelism internally.
 *
 *
 * @section Usage
 *
 * First ectrans needs to be initialized with a function trans_init().
 * This needs to be done only once in the program. It sets up some
 * global structures independent of any resolution.
 *
 * A number of resolutions can now be setup using trans_setup() for each
 * resolution.
 * Every call to trans_setup() involves allocating and computing the
 * transformation coefficients, and should be done only once for
 * every intended resolution as it can be very expensive and requires
 * to store a lot of memory. The resolution can be referred to with
 * a trans "handle" of the Trans_t type.
 *
 * Using this handle, one can now transform fields. Either many fields
 * can be transformed simultaneously, or the transform functions can
 * be called multiple times to transform any number of fields separately.
 *
 * The function to do a transform from gridpoints to spectral is called trans_dirtrans().
 * The function to do a transform from spectral to gridpoints is called trans_invtrans().
 * The function to do the adjoint of the spectral to gridpoints transform is called trans_invtrans_adj().
 * It also transforms the data from gridpoints to spectral.
 *
 * In case of distrubuted parallelism (MPI), the functions trans_dirtrans(), trans_invtrans(),
 * and trans_invtrans_adj() work on distributed fields.
 * In order to convert to and from a global view of the field
 * (e.g. for reading / writing), one can use the functions trans_distspec(), trans_gathspec(),
 * trans_distgrid(), trans_gathgrid().
 *
 * Every handle needs to be cleaned up when no longer required, to release
 * the memory and coefficients stored internally. This can be done with the
 * function trans_delete().
 *
 * Lastly, transi needs to be finalized with trans_finalize(), which will
 * clean up any remaining internal global structures
 *
 * @author Willem Deconinck
 * @date Jul 2014
 */

/*!
 *  @file  transi.h
 *  @brief C-interface to the IFS trans-library
 *
 *  This file declares the C-API to the IFS trans-library
 *  Definitions of routines are implemented in
 *  trans_module.F90, which redirects function calls
 *  to the IFS TRANS library
 *
 *  @author Willem Deconinck (nawd)
 *  @date Jul 2014
 */

#ifndef ectrans_transi_h
#define ectrans_transi_h

#include <stddef.h> // size_t

typedef int _bool;

#ifdef __cplusplus
extern "C" {
#endif

#include "ectrans/version.h"

#define TRANS_FFT992 1
#define TRANS_FFTW   2

#define TRANS_SUCCESS         0

struct Trans_t;
struct DirTrans_t;
struct InvTrans_t;
struct InvTransAdj_t;
struct DistGrid_t;
struct GathGrid_t;
struct DistSpec_t;
struct GathSpec_t;
struct VorDivToUV_t;
struct SpecNorm_t;


/*!
  @brief Get error message relating to error code
 */
const char* trans_error_msg(int errcode);

/*!
  @brief Set limit on maximum simultaneously allocated transforms

  @note Advanced feature

  Default value is 10
  This function needs to be called before trans_init() or trans_setup(), and
  ONLY if the default value needs to be changed.
 */
int trans_set_handles_limit(int limit);

/*!
  @brief Set radius of planet used in trans

  @note Advanced feature

  Default value of radius is Earth's radius (6371.22e+03)
  This function needs to be called before trans_init() or trans_setup(), and
  ONLY if the default value needs to be changed.
 */
int trans_set_radius(double radius);

/*!
  @brief Set nprtrv for parallel distribution of fields in spectral space

  @note Advanced feature

  Default value of nprtrv is 1, meaning that there is no
  parallel distribution of the same wave number for different fields (or levels)
  This function needs to be called before trans_init() or trans_setup(), and
  ONLY if the default value needs to be changed.
 */
int trans_set_nprtrv(int nprtrv);

/*!
  @brief Use MPI in trans library.

  @note Advanced feature

  By default, MPI is used if MPI was detected during compilation.
  To force not to use MPI, this function may be used.
 */
int trans_use_mpi(_bool);

/*!
  @brief Initialize trans library

  This initializes MPI communication, and allocates resolution-independent
  storage. \n
  If this routine is not called, then it will be called internally
  upon the first call to trans_setup()

  @pre call trans_set_radius() and/or trans_set_nprtrv() if radius or
       nprtrv need to be different from default values
 */
int trans_init(void);

int trans_set_read(struct Trans_t*, const char* filepath);
int trans_set_write(struct Trans_t*, const char* filepath);
int trans_set_cache(struct Trans_t*, const void*, size_t);

/*!
  @brief Setup a new resolution to be used in the trans library

  @param trans   Trans_t struct, that needs to have following variables defined:
                   - ndgl   -- number of lattitudes
                   - nloen  -- number of longitudes for each lattitude
                   - nsmax  -- spectral truncation wave number

  All scalar values in the struct will be filled in.
  Remaining array values will be deallocated and set to null.
  To define array values, make individual calls to trans_inquire()

  <b>Usage:</b>
  @code{.c}
  struct Trans_t trans;
  trans_new(&trans);
  trans.ndgl  = ... ;
  trans.nloen = malloc( sizeof(int)*trans.ndgl );
  ... // Read in or compute nloen values
  trans.nsmax = (2*trans.ndgl-1)/2; // For typical linear grid
  trans_setup(&trans);
  @endcode
  @note If trans_init() was not called beforehand, it will be called
        internally
 */
int trans_setup(struct Trans_t* trans);


/*!
  @brief Inquire the trans library for array values

  @param trans   Trans_t struct which needs to have been setup using trans_setup()
  @param varlist comma-separated string of values to inquire

  The inquired values will be allocated if needed, and filled
  in in the Trans_t struct
 */
int trans_inquire(struct Trans_t* trans, const char* varlist);

/*!
  @brief Direct spectral transform (from grid-point to spectral)

  A DirTrans_t struct, initialised with new_dirtrans(),
  groups all arguments

  @param dirtrans  DirTrans_t struct, containing all arguments.

  <b>Usage:</b>
  - Transform of scalar fields
  @code{.c}
  struct Trans_t trans;
  trans_new(&trans);
  ... // Missing setup of trans
  double* rgp       = malloc( sizeof(double) * nscalar*trans.ngptot );
  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  struct DirTrans_t dirtrans = new_dirtrans(&trans);
    dirtrans.nscalar   = nscalar;   // input
    dirtrans.rgp       = rgp;       // input
    dirtrans.rspscalar = rspscalar; // output
  trans_dirtrans(&dirtrans);
  @endcode
  - Transform of U and V fields to vorticity and divergence
  @code
  struct Trans_t trans;
  trans_new(&trans);
  ... // Missing setup of trans
  double* rgp    = malloc( sizeof(double) * 2*nvordiv*trans.ngptot );
  double* rspvor = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspdiv = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  struct DirTrans_t dirtrans = new_dirtrans(&trans);
    dirtrans.nvordiv   = nvordiv;  // input
    dirtrans.rgp       = rgp;      // input    --    order: U, V
    dirtrans.rspvor    = rspvor;   // output
    dirtrans.rspdiv    = rspdiv;   // output
  trans_dirtrans(&dirtrans);
  @endcode
  - Transform of U and V fields to vorticity and divergence, as well as scalar fields
  @code
  struct Trans_t trans;
  trans_new(&trans);
  ... // Missing setup of trans
  double* rgp       = malloc( sizeof(double) * (2*nvordiv+nscalar)*trans.ngptot );
  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  double* rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  struct DirTrans_t dirtrans = new_dirtrans(&trans);
    dirtrans.nvordiv   = nvordiv;   // input
    dirtrans.nscalar   = nscalar;   // input
    dirtrans.rgp       = rgp;       // input   --   order: U, V, scalars
    dirtrans.rspscalar = rspscalar; // output
    dirtrans.rspvor    = rspvor;    // output
    dirtrans.rspdiv    = rspdiv;    // output
  trans_dirtrans(&dirtrans);
  @endcode

  @note trans_dirtrans works on distributed arrays
 */
int trans_dirtrans(struct DirTrans_t* dirtrans);

/*!
  @brief Inverse spectral transform (from spectral grid-point)

  A InvTrans_t struct, initialised with new_invtrans(),
  groups all arguments

  @param invtrans  InvTrans_t struct, containing all arguments.

  <b>Usage:</b>
  - Transform of scalars
  @code{.c}
  struct Trans_t trans;
  trans_new(&trans);
  ... // Missing setup of trans
  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  double* rgp       = malloc( sizeof(double) * nscalar*trans.ngptot );

  // Inverse Transform
  struct InvTrans_t invtrans = new_invtrans(&trans);
    invtrans.nscalar   = nscalar;    // input
    invtrans.rspscalar = rspscalar;  // input
    invtrans.rgp       = rgp;        // output
  trans_invtrans(&invtrans);
  @endcode

  - Transform vorticity and divergence to U and V
  @code{.c}
  struct Trans_t trans;
  trans_new(&trans);
  ... // Missing setup of trans
  double* rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rgp       = malloc( sizeof(double) * 2*nvordiv*trans.ngptot );

  // Inverse Transform
  struct InvTrans_t invtrans = new_invtrans(&trans);
    invtrans.nvordiv   = nvordiv;    // input
    invtrans.rspvor    = rspvor;     // input
    invtrans.rspdiv    = rspdiv;     // input
    invtrans.rgp       = rgp;        // output  --  order: u, v
  trans_invtrans(&invtrans);
  @endcode
  - Transform of vorticity, divergence *and* scalars to U, V, scalars
  @code{.c}
  struct Trans_t trans;
  trans_new(&trans);
  ... // Missing setup of trans
  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  double* rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rgp       = malloc( sizeof(double) * (2*nvordiv+nscalar)*trans.ngptot );

  // Inverse Transform
  struct InvTrans_t invtrans = new_invtrans(&trans);
    invtrans.nscalar   = nscalar;    // input
    invtrans.nvordiv   = nvordiv;    // input
    invtrans.rspscalar = rspscalar;  // input
    invtrans.rspvor    = rspvor;     // input
    invtrans.rspdiv    = rspdiv;     // input
    invtrans.rgp       = rgp;        // output  --  order: u, v, scalars
  trans_invtrans(&invtrans);
  @endcode

  @note trans_invtrans works on distributed arrays
 */
int trans_invtrans(struct InvTrans_t* invtrans);



/*!
  @brief Adjoint of the Inverse spectral transform (from grid-point spectral)

  A InvTransAdj_t struct, initialised with new_invtrans_adj(),
  groups all arguments

  @param invtrans_adj  InvTransAdj_t struct, containing all arguments.

  <b>Usage:</b>
  - Transform of scalars
  @code{.c}
  struct Trans_t trans;
  trans_new(&trans);
  ... // Missing setup of trans
  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  double* rgp       = malloc( sizeof(double) * nscalar*trans.ngptot );

  // Adjoint of Inverse Transform
  struct InvTransAdj_t invtrans_adj = new_invtrans_adj(&trans);
    invtrans_adj.nscalar   = nscalar;    // input
    invtrans_adj.rspscalar = rspscalar;  // output
    invtrans_adj.rgp       = rgp;        // input
  trans_invtrans_adj(&invtrans_adj);
  @endcode

  - Adjoint of Transform vorticity and divergence to U and V
  @code{.c}
  struct Trans_t trans;
  trans_new(&trans);
  ... // Missing setup of trans
  double* rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rgp       = malloc( sizeof(double) * 2*nvordiv*trans.ngptot );

  // Adjoint of Inverse Transform
  struct InvTransAdj_t invtrans_adj = new_invtrans_adj(&trans);
    invtrans_adj.nvordiv   = nvordiv;    // input
    invtrans_adj.rspvor    = rspvor;     // output
    invtrans_adj.rspdiv    = rspdiv;     // output
    invtrans_adj.rgp       = rgp;        // input  --  order: u, v
  trans_invtrans_adj(&invtrans_adj);
  @endcode
  - Adjoint of Transform of vorticity, divergence *and* scalars to U, V, scalars
  @code{.c}
  struct Trans_t trans;
  trans_new(&trans);
  ... // Missing setup of trans
  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  double* rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rgp       = malloc( sizeof(double) * (2*nvordiv+nscalar)*trans.ngptot );

  // Adjoint of Inverse Transform
  struct InvTransAdj_t invtrans_adj = new_invtrans_adj(&trans);
    invtrans_adj.nscalar   = nscalar;    // input
    invtrans_adj.nvordiv   = nvordiv;    // input
    invtrans_adj.rspscalar = rspscalar;  // input
    invtrans_adj.rspvor    = rspvor;     // input
    invtrans_adj.rspdiv    = rspdiv;     // input
    invtrans_adj.rgp       = rgp;        // output  --  order: u, v, scalars
  trans_invtrans_adj(&invtrans_adj);
  @endcode

  @note trans_invtrans_adj works on distributed arrays
 */
int trans_invtrans_adj(struct InvTransAdj_t* invtrans_adj);

/*!
  @brief Distribute global gridpoint array among processors

  <b>Usage:</b>
  @code{.c}
  struct Trans_t trans;
  trans_new(&trans);
  ... // missing setup
  int nfld = 1;
  double* rgpg = NULL;
  if( trans.myproc == 1 ) // Load global field in proc 1
  {
    rgpg = malloc( sizeof(double) * trans.ngptotg*nfld );
    ... // load data in rgpg[nfld][ngptotg]
  }
  int* nfrom = malloc( sizeof(int) * nfld );
  nfrom[0] = 1; // Global field 0 sits in proc 1

  double* rgp  = malloc( sizeof(double) * nfld*trans.ngptot  );
  struct DistGrid_t distgrid = new_distgrid(&trans);
    distgrid.nfrom = nfrom;
    distgrid.rgpg  = rgpg;
    distgrid.rgp   = rgp;
    distgrid.nfld  = nfld;
  trans_distgrid(&distgrid);
  @endcode
 */
int trans_distgrid(struct DistGrid_t* distgrid);

/*!
  @brief Gather global gridpoint array from processors

  <b>Usage:</b>
  @code{.c}
  struct Trans_t trans;
  trans_new(&trans);
  ... // missing setup

  // Distributed field
  int nfld = 1;
  double* rgp  = malloc( sizeof(double) * nfld*trans.ngptot  );
  ... // load data in rgp[nfld][ngptot]

  // Global field
  double* rgpg = NULL;
  if( trans.myproc == 1 ) // We will gather to proc 1
  {
    rgpg = malloc( sizeof(double) * nfld*trans.ngptotg );
  }
  int* nto = malloc( sizeof(int) * nfld );
  nto[0] = 1;

  // Gather global fields
  struct GathGrid_t gathgrid = new_gathgrid(&trans);
    gathgrid.rgp  = rgp;
    gathgrid.rgpg = rgpg;
    gathgrid.nto  = nto;
    gathgrid.nfld = nfld;
  trans_gathgrid(&gathgrid);
  @endcode
 */
int trans_gathgrid(struct GathGrid_t* gathgrid);

/*!
  @brief Distribute global spectral array among processors

  <b>Usage:</b>
  @code{.c}
  struct Trans_t trans;
  trans_new(&trans);
  ... // missing setup

  // Global fields to be distributed
  int nscalar = 1;
  double* rspscalarg = NULL;
  if( trans.myproc == 1 )
  {
    rspscalarg = malloc( sizeof(double) * nscalar*trans.nspec2g );
    ... // load data in rspscalarg[nspec2g][nscalar]
  }
  int* nfrom = malloc( sizeof(int) * nscalar );
  nfrom[0] = 1;

  // Distribute to local fields
  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  struct DistSpec_t distspec = new_distspec(&trans);
    distspec.rspec  = rspscalar;
    distspec.rspecg = rspscalarg;
    distspec.nfld   = nscalar;
    distspec.nfrom  = nto;
  trans_distspec(&distspec);
  @endcode
 */
int trans_distspec(struct DistSpec_t* distspec);

/*!
   @brief Gather global spectral array from processors

   <b>Usage:</b>
   @code{.c}
   struct Trans_t trans;
   trans_new(&trans);
   ... // missing setup

   // We have some distributed spectral fields "rspscalar"
   int nscalar = 1;
   double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
   ... // load data in rspscalar[nspec2][nscalar]

   // We want to gather to proc 1
   double* rspscalarg = NULL;
   if( trans.myproc == 1 )
     rspscalarg = malloc( sizeof(double) * nscalar*trans.nspec2g );
   int* nto = malloc( sizeof(int) * nscalar );
   nto[0] = 1;
   struct GathSpec_t gathspec = new_gathspec(&trans);
     gathspec.rspec  = rspscalar;
     gathspec.rspecg = rspscalarg;
     gathspec.nfld   = nscalar;
     gathspec.nto    = nto;
   trans_gathspec(&gathspec);
   @endcode
 */
int trans_gathspec(struct GathSpec_t* gathspec);

/*!
   @brief Convert Spectral vorticity & divergence to Spectral u*cos(theta) & v*cos(theta)

   <b>Usage:</b>
   @code{.c}
   // We have some global spectral fields for vorticity,divergence,u*cos(theta),v*cos(theta)
   int nfld = 1;
   double* rspvor = malloc( sizeof(double) * nfld*ncoeff );
   double* rspdiv = malloc( sizeof(double) * nfld*ncoeff );
   double* rspu   = malloc( sizeof(double) * nfld*ncoeff );
   double* rspv   = malloc( sizeof(double) * nfld*ncoeff );
   ... // load data in rspvor[ncoeff][nfld]
   ... // load data in rspdiv[ncoeff][nfld]

   struct VorDivToUV_t vordiv_to_UV = new_vordiv_to_UV();
     vordiv_to_UV.rspvor = rspvor;
     vordiv_to_UV.rspdiv = rspdiv;
     vordiv_to_UV.rspu   = rspu;
     vordiv_to_UV rspv   = rspv;
     vordiv_to_UV.nfld   = nfld;
     vordiv_to_UV.ncoeff = ncoeff;
     vordiv_to_UV.nsmax  = nsmax;
   trans_vordiv_to_UV(&vordiv_to_UV);
   @endcode

   @note
   - nfld indicates the multiplicity for each variable seperately
   - ncoeff is equivalent to trans.nspec2 for distributed, and trans.nspec2g for global fields
   - nsmax indicates the spectral truncation T.
 */
int trans_vordiv_to_UV(struct VorDivToUV_t* vordiv_to_UV);

/*!
  @brief Compute global spectral norms

  <b>Usage:<\b>
  @code{.c}
  int nfld = 1;
  double* rspec = malloc( sizeof(double) * nfld*trans.nspec2 );
  double* rnorm = malloc( sizeof(double) * nfld );
  ... // load data in rspec[nspec2][nfld]

  struct SpecNorm_t specnorm = new_specnorm(&trans);
    specnorm.rspec = rspec;
    specnorm.rnorm = rnorm;
    specnorm.nfld  = nfld;
  trans_specnorm(specnorm);
  @endcode
*/
int trans_specnorm(struct SpecNorm_t* specnorm);

/*!
   @brief Remove footprint of specific resolution

   @param trans  Trans_t struct describing specific resolution

   All arrays will be deallocated.
 */
int trans_delete(struct Trans_t* trans);


/*!
  @brief Finalize trans library

  This finalizes MPI communication, and deallocates resolution-independent
  storage. After this, no more calls to trans should be made
 */
int trans_finalize(void);


/*!
   @brief Struct that holds information to do transforms
          for one particular grid resolution

   The values ndgl, nloen, and nsmax need to be provided yourself, all other
   values will be defined during the trans_setup() call or trans_inquire() calls

   - All scalar values will be defined by trans_setup()
   - All array values will be allocated if needed, and defined by
     individual calls to trans_inquire()

   @note Many of these values are of no interest for normal usage
 */
struct Trans_t {

  /*! @{ @name INPUT */
  int    ndgl;       //!< @brief  Number of lattitudes
  int*   nloen;      //!< @brief  Number of longitude points for each latitude \n
                     //!<         DIMENSIONS(1:NDGL)
  int    nlon;       //!< @brief  Number of longitude points for all latitudes \n
  int    nsmax;      //!< @brief  Spectral truncation wave number

  _bool  lsplit;         //!< @brief  If false, the distribution does not allow latitudes to be split
  int    llatlon;        //!< @brief  If true, the transforms compute extra coefficients for
                         //!<         latlon transforms
  int    flt;            //!< @brief  If true, the Fast-Legendre-Transform method is used
                         //!<         which is faster for higher resolutions (N1024)
  int    fft;            //!< @brief  FFT library to use underneith \n
                         //!<         FFT992 = 1 ; FFTW = 2

  char*  readfp;
  char*  writefp;
  const void*  cache;
  size_t cachesize;
  /*! @} */

  /*! @{ @name PARALLELISATION */
  int    myproc;     //!< @brief  Current MPI task (numbering starting at 1)
  int    nproc;      //!< @brief  Number of parallel MPI tasks
  /*! @} */

  /*! @{ @name MULTI-TRANSFORMS-MANAGEMENT */
  int   handle;      //!< @brief  Resolution tag for which info is required ,default is the
                     //!<         first defined resulution (input)
  _bool ldlam;       //!< @brief  True if the corresponding resolution is LAM, false if it is global
  /*! @} */

  /*! @{ @name SPECTRAL SPACE */
  int   nspec;       //!< @brief  Number of complex spectral coefficients on this PE
  int   nspec2;      //!< @brief  Number of complex spectral coefficients on this PE times 2 (real and imag)
  int   nspec2g;     //!< @brief  global KSPEC2
  int   nspec2mx;    //!< @brief  Maximun KSPEC2 among all PEs
  int   nump;        //!< @brief  Number of spectral waves handled by this PE
  int   ngptot;      //!< @brief  Total number of grid columns on this PE
  int   ngptotg;     //!< @brief  Total number of grid columns on the Globe
  int   ngptotmx;    //!< @brief  Maximum number of grid columns on any of the PEs
  int*  ngptotl;     //!< @brief  Number of grid columns on each PE \n
                     //!<         DIMENSIONS(1:N_REGIONS_NS,1:N_REGIONS_EW)
  int*  nmyms;       //!< @brief  This PEs spectral zonal wavenumbers
                     //!<         DIMENSIONS(1:NUMP)
  int*  nasm0;       //!< @brief  Address in a spectral array of (m, n=m) \n
                     //!<         DIMENSIONS(0:NSMAX)
  int   nprtrw;      //!< @brief  Number of processors in A-direction (input)
  int*  numpp;       //!< @brief  No. of wave numbers each wave set is responsible for. \n
                     //!<         DIMENSIONS(1:NPRTRW)
  int*  npossp;      //!< @brief  Defines partitioning of global spectral fields among PEs \n
                     //!<         DIMENSIONS(1:NPRTRW+1)
  int*  nptrms;      //!< @brief  Pointer to the first wave number of a given a-set \n
                     //!<         DIMENSIONS(1:NPRTRW)
  int*  nallms;      //!< @brief  Wave numbers for all wave-set concatenated together
                     //!<         to give all wave numbers in wave-set order \n
                     //!<         DIMENSIONS(1:NSMAX+1)
  int*  ndim0g;      //!< @brief  Defines partitioning of global spectral fields among PEs \n
                     //!<         DIMENSIONS(0:NSMAX)
  int*  nvalue;      //!< @brief  n value for each KSPEC2 spectral coeffient\n
                     //!<         DIMENSIONS(1:NSPEC2)
  /*! @} */

  /*! @{ @name GRIDPOINT SPACE */
  int   n_regions_NS;//!< @brief   Number of regions in North-South direction
  int   n_regions_EW;//!< @brief   Number of regions in East-West direction
  int   my_region_NS;//!< @brief   My region in North-South direction
  int   my_region_EW;//!< @brief   My region in East-West direction
  int*  n_regions;   //!< @brief   Number of East-West Regions per band of North-South Regions
                     //!< @brief   DIMENSIONS(1:N_REGIONS_NS)
  int*  nfrstlat;    //!< @brief   First latitude of each a-set in grid-point space
                     //!<          DIMENSIONS(1:N_REGIONS_NS)
  int*  nlstlat;     //!< @brief   Last latitude of each a-set in grid-point space
                     //!<          DIMENSIONS(1:N_REGIONS_NS)
  int   nfrstloff;   //!< @brief   Offset for first lat of own a-set in grid-point space
  int*  nptrlat;     //!< @brief   Pointer to the start of each latitude
                     //!<          DIMENSIONS(1:NDGL)
  int*  nptrfrstlat; //!< @brief   Pointer to the first latitude of each a-set in
                     //!<          NSTA and NONL arrays
                     //!<          DIMENSIONS(1:N_REGIONS_NS)
  int*  nptrlstlat;  //!< @brief   Pointer to the last latitude of each a-set in
                     //!<          NSTA and NONL arrays
                     //!<          DIMENSIONS(1:N_REGIONS_NS)
  int   nptrfloff;   //!< @brief   Offset for pointer to the first latitude of own a-set
                     //!<          NSTA and NONL arrays, i.e. nptrfrstlat(myseta)-1
  int*  nsta;        //!< @brief   Position of first grid column for the latitudes on a
                     //!<          processor. \n
                     //!<          DIMENSIONS(1:NDGL+N_REGIONS_NS-1,1:N_REGIONS_EW)
                     //!<          @details The information is available for all processors.
                     //!<          The b-sets are distinguished by the last dimension of
                     //!<          nsta(). The latitude band for each a-set is addressed by
                     //!<          nptrfrstlat(jaset),nptrlstlat(jaset), and
                     //!<          nptrfloff=nptrfrstlat(myseta) on this processors a-set.
                     //!<          Each split latitude has two entries in nsta(,:) which
                     //!<          necessitates the rather complex addressing of nsta(,:)
                     //!<          and the overdimensioning of nsta by N_REGIONS_NS.
  int*  nonl;        //!< @brief   Number of grid columns for the latitudes on a processor.
                     //!<          Similar to nsta() in data structure. \n
                     //!<          DIMENSIONS(1:NDGL+N_REGIONS_NS-1,1:N_REGIONS_EW)
  _bool* ldsplitlat; //!< @brief   True if latitude is split in grid point space over
                     //!<          two a-sets. \n
                     //!<          DIMENSIONS(1:NDGL)
  /*! @} */

  /*! @{ @name FOURIER SPACE */
  int   nprtrns;     //!< @brief  No. of sets in N-S direction (Fourier space)
                     //!<         (always equal to NPRTRW)
  int*  nultpp;      //!< @brief  Number of latitudes for which each a-set is calculating
                     //!<         the FFT's. \n
                     //!<         DIMENSIONS(1:NPRTRNS)
  int*  nptrls;      //!< @brief  Pointer to first global latitude of each a-set for which
                     //!<         it performs the Fourier calculations \n
                     //!<         DIMENSIONS(1:NPRTRNS)
  int*  nnmeng;      //!< @brief  associated (with NLOENG) cut-off zonal wavenumber \n
                     //!<         DIMENSIONS(1:NDGL)
  /*! @} */


  /*! @{ @name LEGENDRE  */
  double*  rmu;      //!< @brief  sin(Gaussian latitudes) \n
                     //!<         DIMENSIONS(1:NDGL)
  double*  rgw;      //!< @brief  Gaussian weights \n
                     //!<         DIMENSIONS(1:NDGL)
  double*  rpnm;     //!< @brief  Legendre polynomials \n
                     //!<         DIMENSIONS(1:NLEI3,1:NSPOLEGL)
  int      nlei3;    //!< @brief  First dimension of Legendre polynomials
  int      nspolegl; //!< @brief  Second dimension of Legendre polynomials
  int*     npms;     //!< @brief  Adress for legendre polynomial for given M (NSMAX) \n
                     //!<         DIMENSIONS(0:NSMAX)
  double*  rlapin;   //!< @brief  Eigen-values of the inverse Laplace operator \n
                     //!<         DIMENSIONS(-1:NSMAX+2)
  int*     ndglu;    //!< @brief  Number of active points in an hemisphere for a given wavenumber "m" \n
                     //!<         DIMENSIONS(0:NSMAX)
  /*! @} */

};

/*!
   @brief Constructor for Trans_t, setting default values
   @return Trans_t struct to be used as argument for trans_setup()
 */
int trans_new( struct Trans_t* );

/*!
   @brief Set gridpoint resolution for trans
   @param  trans [in] Trans_t used to setup
   @param  ndgl  [in] Number of lattitudes
   @param  nloen [in] Number of longitude points for each latitude \n
                      DIMENSIONS(1:NDGL)
 */
int trans_set_resol( struct Trans_t* trans, int ndgl, const int* nloen );

/*!
   @brief Set gridpoint resolution for trans for longitude-latitude grids
   @param  trans [in] Trans_t used to setup
   @param  nlon  [in] Number of longitudes
   @param  nlat  [in] Number of latitudes (pole to pole)

   - If nlat is odd, the grid must includes poles and equator
   - If nlat is even, the grid must be its dual (excluding pole and equator),
       so points are shifted with 0.5*dx and 0.5*dy
 */
int trans_set_resol_lonlat( struct Trans_t* trans, int nlon, int nlat );

/*!
   @brief Set spectral truncation wave number for trans
   @param  trans [in] Trans_t used to setup
   @param  nsmax [in] Spectral truncation wave number
 */
int trans_set_trunc( struct Trans_t* trans, int nsmax );


/*!
   @brief Arguments structure for trans_dirtrans()

   Use new_dirtrans() to initialise defaults for the struct (constructor)
 */
struct DirTrans_t
{
  const double* rgp;     //!< @brief  [input] gridpoint fields
                         //!<         @details Dimensioning:  rgp[#ngpblks][2*#nvordiv+#nscalar][#nproma]\n\n
                         //!<         The ordering of the output fields is as follows (all
                         //!<         parts are optional depending on the input switches):
                         //!<         - u       : if #nvordiv > 0
                         //!<         - v       : if #nvordiv > 0
                         //!<         - scalars : if #nscalar > 0
  double* rspscalar;     //!< @brief  [output] spectral scalar valued fields
                         //!<         @details Dimensioning: rspscalar[@link Trans_t::nspec2 nspec2 @endlink][#nscalar]
  double* rspvor;        //!< @brief  [output] spectral vorticity
                         //!<         @details Dimensioning: rspvor[@link Trans_t::nspec2 nspec2 @endlink][#nvordiv]
  double* rspdiv;        //!< @brief  [output] spectral divergence
                         //!<         @details Dimensioning: rspvor[@link Trans_t::nspec2 nspec2 @endlink][#nvordiv]
  int nproma;            //!< @brief  [input,default=@link Trans_t::ngptot ngptot@endlink] Blocking factor for distributed gridpoint array
  int nscalar;           //!< @brief  [input,default=0] Number of scalar fields present in RGP
  int nvordiv;           //!< @brief  [input,default=0] Number of vorticity/divergence fields in RGP
  int ngpblks;           //!< @brief  [input,default=1] Blocking factor for distributed gridpoint array
  int lglobal;           //!< @brief  [input,default=0] rgp is a global input field --> nproma==1,ngpblks==ngptotg
  struct Trans_t* trans; //!< @brief Internal storage of trans object
  int count;             //!< @brief Internal storage for calls to trans_dirtrans()
};
/*!
   @brief Constructor for DirTrans_t, resetting default values
   @param  trans [in] Trans_t used to set defaults
   @return DirTrans_t struct to be used as argument for trans_dirtrans()
 */
struct DirTrans_t new_dirtrans(struct Trans_t* trans);

/*!
   @brief Arguments structure for trans_invtrans()

   Use new_invtrans() to initialise defaults for the struct (constructor)
 */
struct InvTrans_t
{
  const double* rspscalar;  //!< @brief  [input,default=NULL] spectral scalar valued fields
                            //!<         @details Dimensioning: rspscalar[@link Trans_t::nspec2 nspec2 @endlink][#nscalar]
  const double* rspvor;     //!< @brief  [input] spectral vorticity
                            //!<         @details Dimensioning: rspvor[@link Trans_t::nspec2 nspec2 @endlink][#nvordiv]
  const double* rspdiv;     //!< @brief  [input] spectral divergence
                            //!<         @details Dimensioning: rspvor[@link Trans_t::nspec2 nspec2 @endlink][#nvordiv]
  double* rgp;              //!< @brief  [output] gridpoint fields
                            //!<         @details Dimensioning:  rgp[#ngpblks][2*#nvordiv+#nscalar][#nproma]\n\n
                            //!<         The ordering of the output fields is as follows (all
                            //!<         parts are optional depending on the input switches):
                            //!          - vorticity  : if #nvordiv > 0  and  #lvordivgp true
                            //!          - divergence : if #nvordiv > 0  and  #lvordivgp true
                            //!<         - u          : if #nvordiv > 0
                            //!<         - v          : if #nvordiv > 0
                            //!<         - scalars    : if #nscalar > 0
                            //!<         - N-S derivative of scalars : if #nscalar > 0  and  #lscalarders true
                            //!<         - E-W derivative of u       : if #nvordiv > 0  and  #luvders true
                            //!<         - E-W derivative of v       : if #nvordiv > 0  and  #luvders true
                            //!<         - E-W derivative of scalars : if #nscalar > 0  and  #lscalarders true
  int nproma;               //!< @brief  [input,default=@link Trans_t::ngptot ngptot@endlink] Blocking factor for distributed gridpoint array
  int nscalar;              //!< @brief  [input,default=0] Number of scalar fields present in RGP
  int nvordiv;              //!< @brief  [input,default=0] Number of vorticity/divergence fields in RGP
  int lscalarders;          //!< @brief  [input,default=0] Indicate if derivatives of scalars are requested
  int luvder_EW;            //!< @brief  [input,default=0] Indicate if East-West derivative of u and v is requested
  int lvordivgp;            //!< @brief  [input,default=0] Indicate if grid-point vorticity and divergence is requested
  int ngpblks;              //!< @brief  [input,default=1] Blocking factor for distributed gridpoint array
  int lglobal;              //!< @brief  [input,default=0] rgp is a global output field --> nproma==1,ngpblks==ngptotg
  struct Trans_t* trans;    //!< @brief Internal storage of trans object
  int count;                //!< @brief Internal storage for calls to trans_invtrans()
};
/*!
 * @brief Constructor for InvTrans_t, resetting default values
 * @param  trans [in] Trans_t used to set defaults
 * @return InvTrans_t struct to be used as argument for trans_invtrans()
 */
struct InvTrans_t new_invtrans(struct Trans_t* trans);


//! Adjoint of spectral inverse.

struct InvTransAdj_t
{
  double* rspscalar;        //!< @brief  [output,default=NULL] spectral scalar valued fields
                            //!<         @details Dimensioning: rspscalar[@link Trans_t::nspec2 nspec2 @endlink][#nscalar]
  double* rspvor;           //!< @brief  [output] spectral vorticity
                            //!<         @details Dimensioning: rspvor[@link Trans_t::nspec2 nspec2 @endlink][#nvordiv]
  double* rspdiv;           //!< @brief  [output] spectral divergence
                            //!<         @details Dimensioning: rspvor[@link Trans_t::nspec2 nspec2 @endlink][#nvordiv]
  const double* rgp;        //!< @brief  [input] gridpoint fields
                            //!<         @details Dimensioning:  rgp[#ngpblks][2*#nvordiv+#nscalar][#nproma]\n\n
                            //!<         The ordering of the output fields is as follows (all
                            //!<         parts are optional depending on the input switches):
                            //!          - vorticity  : if #nvordiv > 0  and  #lvordivgp true
                            //!          - divergence : if #nvordiv > 0  and  #lvordivgp true
                            //!<         - u          : if #nvordiv > 0
                            //!<         - v          : if #nvordiv > 0
                            //!<         - scalars    : if #nscalar > 0
                            //!<         - N-S derivative of scalars : if #nscalar > 0  and  #lscalarders true
                            //!<         - E-W derivative of u       : if #nvordiv > 0  and  #luvders true
                            //!<         - E-W derivative of v       : if #nvordiv > 0  and  #luvders true
                            //!<         - E-W derivative of scalars : if #nscalar > 0  and  #lscalarders true
  int nproma;               //!< @brief  [input,default=@link Trans_t::ngptot ngptot@endlink] Blocking factor for distributed gridpoint array
  int nscalar;              //!< @brief  [input,default=0] Number of scalar fields present in RGP
  int nvordiv;              //!< @brief  [input,default=0] Number of vorticity/divergence fields in RGP
  int lscalarders;          //!< @brief  [input,default=0] Indicate if derivatives of scalars are requested
  int luvder_EW;            //!< @brief  [input,default=0] Indicate if East-West derivative of u and v is requested
  int lvordivgp;            //!< @brief  [input,default=0] Indicate if grid-point vorticity and divergence is requested
  int ngpblks;              //!< @brief  [input,default=1] Blocking factor for distributed gridpoint array
  int lglobal;              //!< @brief  [input,default=0] rgp is a global output field --> nproma==1,ngpblks==ngptotg
  struct Trans_t* trans;    //!< @brief Internal storage of trans object
  int count;                //!< @brief Internal storage for calls to trans_invtrans_adj()
};
/*!
 * @brief Constructor for InvTransAdj_t, resetting default values
 * @param  trans [in] Trans_t used to set defaults
 * @return InvTransAdj_t struct to be used as argument for trans_invtrans_adj()
 */
struct InvTransAdj_t new_invtrans_adj(struct Trans_t* trans);



/*!
   @brief Arguments structure for trans_distgrid()

   Use new_distgrid() to initialise defaults for the struct (constructor)
 */
struct DistGrid_t
{
  const double* rgpg;    //!< @brief  Global gridpoint array
                         //!<         Fortran DIMENSIONS(1:NGPTOTG,1:NFLDG)
                         //!<         C/C++   DIMENSIONS[NFLDG][NGPTOTG]
  double* rgp;           //!< @brief  Distributed gridpoint array
                         //!<         Fortran DIMENSIONS(1:NPROMA,1:NFLD,1:NGPBLKS)
                         //!<         C/C++   DIMENSIONS[NGPBLKS][NFDL][NPROMA]
                         //!<         Default: NPROMA=NGPTOT, NGPBLKS=1
  const int* nfrom;      //!< @brief  Processors responsible for distributing each field
                         //!<         DIMENSIONS(1:NFLD)
  int nproma;            //!< @brief  Blocking factor for distributed gridpoint array
  int nfld;              //!< @brief  Number of distributed fields
  int ngpblks;           //!< @brief  Blocking factor for distributed gridpoint array
  struct Trans_t* trans; //!< @brief Internal storage of trans object
  int count;             //!< @brief Internal storage for calls to trans_invtrans()
};
/*!
   @brief Constructor for DistGrid_t, resetting default values
   @param  trans [in] Trans_t used to set defaults
   @return DistGrid_t struct to be used as argument for trans_distgrid()
 */
struct DistGrid_t new_distgrid(struct Trans_t* trans);

/*!
   @brief Arguments structure for trans_gathgrid()

   Use new_gathgrid() to initialise defaults for the struct (constructor)
 */
struct GathGrid_t
{
  double* rgpg;          //!< @brief  Global gridpoint array
                         //!<         Fortran DIMENSIONS(1:NGPTOTG,1:NFLDG)
                         //!<         C/C++   DIMENSIONS[NFLDG][NGPTOTG]
                         //!<         DIMENSIONS(1:NFLDG,1:NGPTOTG)
  const double* rgp;     //!< @brief  Distributed gridpoint array
                         //!<         Fortran DIMENSIONS(1:NPROMA,1:NFLD,1:NGPBLKS)
                         //!<         C/C++   DIMENSIONS[NGPBLKS][NFDL][NPROMA]
                         //!<         Default: NPROMA=NGPTOT, NGPBLKS=1
  const int* nto;        //!< @brief  Processors responsible for gathering each field
                         //!<         Fortran DIMENSIONS(1:NFLD)
  int nproma;            //!< @brief  Blocking factor for distributed gridpoint array
  int nfld;              //!< @brief  Number of distributed fields
  int ngpblks;           //!< @brief  Blocking factor for distributed gridpoint array
  struct Trans_t* trans; //!< @brief Internal storage of trans object
  int count;             //!< @brief Internal storage for calls to trans_invtrans()
};
/*!
  @brief Constructor for GathGrid_t, resetting default values
  @param  trans [in] Trans_t used to set defaults
  @return GathGrid_t struct to be used as argument for trans_gathgrid()
 */
struct GathGrid_t new_gathgrid(struct Trans_t* trans);

/*!
   @brief Arguments structure for trans_distspec()

   Use new_distspec() to initialise defaults for the struct (constructor)
 */
struct DistSpec_t
{
  const double* rspecg;  //!< @brief  Global spectral array
                         //!<         Fortran DIMENSIONS(1:NFLDG,1:NSPEC2G)
                         //!<         C/C++   DIMENSIONS[NSPEC2G][NFLDG]
  double* rspec;         //!< @brief  Local spectral array
                         //!<         Fortran DIMENSIONS(1:NFLD,1:NSPEC2)
                         //!<         C/C++   DIMENSIONS[NSPEC2][NFLD]
  const int* nfrom;      //!< @brief  Processors responsible for distributing each field
                         //!<         Fortran DIMENSIONS(1:NFLD)
  int nfld;              //!< @brief  Number of distributed fields
  struct Trans_t* trans; //!< @brief Internal storage of trans object
  int count;             //!< @brief Internal storage for calls to trans_invtrans()
};
/*!
   @brief Constructor for DistSpec_t, resetting default values
   @param  trans [in] Trans_t used to set defaults
   @return DistSpec_t struct to be used as argument for trans_distspec()
 */
struct DistSpec_t new_distspec(struct Trans_t* trans);

/*!
   @brief Arguments structure for trans_gathspec()

   Use new_gathspec() to initialise defaults for the struct (constructor)
 */
struct GathSpec_t
{
  double* rspecg;        //!< @brief  Global spectral array
                         //!<         Fortran DIMENSIONS(1:NFLDG,1:NSPEC2G)
                         //!<         C/C++   DIMENSIONS[NSPEC2G][NFLDG]
  const double* rspec;   //!< @brief  Local spectral array
                         //!<         Fortran DIMENSIONS(1:NFLD,1:NSPEC2)
                         //!<         C/C++   DIMENSIONS[NSPEC2][NFLD]
  const int* nto;        //!< @brief  Processors responsible for gathering each field
                         //!<         DIMENSIONS(1:NFLD)
  int nfld;              //!< @brief  Number of distributed fields
  struct Trans_t* trans; //!< @brief Internal storage of trans object
  int count;             //!< @brief Internal storage for calls to trans_invtrans()
};
/*!
   @brief Constructor for GathSpec_t, resetting default values
   @param  trans [in] Trans_t used to set defaults
   @return GathSpec_t struct to be used as argument for trans_gathspec()
 */
struct GathSpec_t new_gathspec(struct Trans_t* trans);

/*!
   @brief Arguments structure for trans_vordiv_to_UV()

   Use new_vordiv_to_uv() to initialise defaults for the struct (constructor)
 */
struct VorDivToUV_t
{
  const double* rspvor;  //!< @brief  Local spectral array for vorticity
                         //!<         Fortran DIMENSIONS(1:NFLD,1:NSPEC2)
                         //!<         C/C++   DIMENSIONS[NSPEC2][NFLD]
  const double* rspdiv;  //!< @brief  Local spectral array for divergence
                         //!<         Fortran DIMENSIONS(1:NFLD,1:NSPEC2)
                         //!<         C/C++   DIMENSIONS[NSPEC2][NFLD]
  double* rspu;          //!< @brief  Local spectral array for U*cos(theta)
                         //!<         Fortran DIMENSIONS(1:NFLD,1:NSPEC2)
                         //!<         C/C++   DIMENSIONS[NSPEC2][NFLD]
  double* rspv;          //!< @brief  Local spectral array for V*cos(theta)
                         //!<         Fortran DIMENSIONS(1:NFLD,1:NSPEC2)
                         //!<         C/C++   DIMENSIONS[NSPEC2][NFLD]
  int nfld;              //!< @brief  Number of distributed fields
  int nsmax;             //!< @brief  Spectral resolution (T)
  int ncoeff;            //!< @brief  number of spectral coefficients
                         //!<         (equivalent to nspec2 for distributed or nspec2g for global)
  int count;             //!< @brief Internal storage for calls to trans_vordiv_toUV()
};
/*!
   @brief Constructor for VorDivToUV_t, resetting default values
   @return VorDivToUV_t struct to be used as argument for trans_gathspec()
 */
struct VorDivToUV_t new_vordiv_to_UV(void);


struct SpecNorm_t
{
  const double *rspec;     //!< @brief Spectral array to compute norm of
                           //!<         Fortran DIMENSIONS(1:NFLD,1:NSPEC2)
                           //!<         C/C++   DIMENSIONS[NSPEC2][NFLD]
  int nmaster;             //!< @brief Processor to receive norms (value 1 means MPI_RANK 0)
  const double *rmet;      //!< @brief metric, OPTIONAL
                           //!         DIMENSIONS(0:NSMAX)
  double* rnorm;           //!< @brief Norms (output for processor nmaster)
                           //!<        DIMENSIONS(1:NFLD)
  int nfld;                //!< @brief Number of fields
  struct Trans_t* trans;   //!< @brief Internal storage of trans object
  int count;               //!< @brief Internal storage for calls to trans_invtrans()
};
/*!
   @brief Constructor for SpecNorm_t, resetting default values
   @return SpecNorm_t struct to be used as argument for trans_specnorm()
 */
struct SpecNorm_t new_specnorm(struct Trans_t* trans);

#ifdef __cplusplus
}
#endif

#endif
