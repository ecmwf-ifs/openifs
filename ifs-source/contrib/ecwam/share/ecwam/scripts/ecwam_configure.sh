# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

function ecwam_configure {
  #wave_set_config version:  20120305
  set -eu
  set +x
  ########################################################################
  #                                                                      #
  #     function ecwam_configure                                         #
  #     ------------------------                                         #
  #                                                                      #
  #     Purpose     : sets variable for the wave model configurations:   #
  #     -------                                                          # 
  #     xdella : grid spacing in degree                                  #
  #     amosop : most southern boundary                                  #
  #     amonop : most northern boundary                                  #
  #     amowep : most western boundary                                   # 
  #     amoeap : most eastern boundary                                   #
  #     iper : iper=1 for periodic domain else = 0                       #
  #     irgg : irgg=1 for irregular lat-lon grid else = 0 for regular    #
  #     llobstrct: T if the bahtymetry and the unresolved bathymetry     #
  #                parametrisation are computed, F if old wam type of    #
  #                mean bathymetry data are used instead.                #
  #     fr1 : value of the reference frequency wrt ifre1 (below)         #
  #     ifre1 : index of the reference frequency in the freqency array   #
  #     deptha: water minimum depth for depth tables (metres)            #
  #     cldomain: g for global model, m for limited area model and       #
  #               s for swamp and one grid point cases.                  #
  #     llunstr   unstructured grid option
  #                                                                      #
  #     Usage   : ecwam_configure $wamresol $wamnfre $wambathy           #
  #     -----                                                            #
  #                                                                      #
  ########################################################################

  wamresol=$1
  wamnfre=$2
  wambathy=$3

  # valid for all configurations, except if stated otherwise

  fr1=4.177248E-02
  if [[ $wamnfre = 25 ]] ; then
    ifre1=1
  else
    ifre1=3
  fi
  if [[ $wambathy = "ETOPO1" ]] ; then
    deptha=3.0
  else
    deptha=5.0
  fi

  cldomain=g
  llobstrct=T

  llunstr=F
  irgg=1

  # set up different configuration
  if [[ ${wamresol} =~ ^O[0-9]+$ ]]; then  # Octahedral grids pattern.   (Grid for TCO1279: O1280)
    ${SCRIPTS_DIR}/ecwam_grids.py ${wamresol} --grid_description > grid_description
    cwamresol=$(${SCRIPTS_DIR}/ecwam_grids.py ${wamresol} --variable=cwamresol)
    xdella=$(${SCRIPTS_DIR}/ecwam_grids.py ${wamresol} --variable=dlat)
    xdello=$(${SCRIPTS_DIR}/ecwam_grids.py ${wamresol} --variable=dlon)
    amonop=$(${SCRIPTS_DIR}/ecwam_grids.py ${wamresol} --variable=north)
    amosop=$(${SCRIPTS_DIR}/ecwam_grids.py ${wamresol} --variable=south)
    amowep=$(${SCRIPTS_DIR}/ecwam_grids.py ${wamresol} --variable=west)
    amoeap=$(${SCRIPTS_DIR}/ecwam_grids.py ${wamresol} --variable=east)
    irgg=$(${SCRIPTS_DIR}/ecwam_grids.py ${wamresol} --variable=irgg)
    iper=1
  elif [[ $wamresol = global5 ]] ; then
    xdella=0.0500000
    amosop=-78.00000
    amonop=90.00000
    amowep=0.0
    amoeap=359.9500
    iper=1
  elif [[ $wamresol = global10 ]] ; then
    xdella=0.1000000
    amosop=-78.00000
    amonop=90.00000
    amowep=0.0
    amoeap=359.9000
    iper=1
  elif [[ $wamresol = global12.5 ]] ; then
    xdella=0.1250000
    amosop=-78.00000
    amonop=90.00000
    amowep=0.0
    amoeap=359.8750
    iper=1
  elif [[ $wamresol = global25 ]] ; then
    xdella=0.2500000
    amosop=-78.00000
    amonop=90.00000
    amowep=0.0
    amoeap=359.7500
    iper=1
  elif [[ $wamresol = global36 ]] ; then
    xdella=0.3600000
    amosop=-78.12000
    amonop=90.00000
    amowep=0.0
    amoeap=359.6400
    iper=1
  elif [[ $wamresol = global50 ]] ; then
    xdella=0.5000000
    amosop=-78.00000
    amonop=90.00000
    amowep=0.0
    amoeap=359.5000
    iper=1
  elif [[ $wamresol = global100 ]] ; then
    xdella=1.0000000
    amosop=-78.00000
    amonop=90.00000
    amowep=0.0
    amoeap=359.0000
    iper=1
  elif [[ $wamresol = global150 ]] ; then
    xdella=1.5000000
    amosop=-78.00000
    amonop=90.00000
    amowep=0.0
    amoeap=358.5000
    deptha=5.0
    iper=1
  elif [[ $wamresol = global300 ]] ; then
    xdella=3.0000000
    amosop=-78.00000
    amonop=90.00000
    amowep=0.0
    amoeap=357.0000
    deptha=5.0
    iper=1
  elif [[ $wamresol = medite25 ]] ; then
    xdella=0.2500000
    amosop=9.00000
    amonop=90.00000
    amowep=-98.00000
    amoeap=42.0000
    iper=0
    ##!!
    cldomain=m
    ##!!
  elif [[ $wamresol = medite25 ]] ; then
    xdella=0.2500000
    amosop=30.00000
    amonop=46.00000
    amowep=-6.00000
    amoeap=36.0000
    iper=0
    fr1=0.05
    ifre1=1
    ##!!
    cldomain=m
    ##!!
  elif [[ $wamresol = medite15 ]] ; then
    xdella=0.1500000
    amosop=5.10000
    amonop=90.00000
    amowep=-98.00000
    amoeap=53.9500
    iper=0
    ##!!
    cldomain=m
    ##!!
  elif [[ $wamresol = medite10 ]] ; then
    xdella=0.1000000
    amosop=-78.00000
    amonop=90.00000
    amowep=0.0
    amoeap=359.9000
    deptha=2.0
    iper=1
    ##!!
    cldomain=m
    ##!!
  elif [[ $wamresol = onegrdpt ]] ; then
  # one grid point setup
    xdella=0.5000000
    amosop=0.00000
    amonop=0.00000
    amowep=2.50000
    amoeap=2.5000
    iper=0
    deptha=1.0
    irgg=0
    ##!!
    cldomain=s
    llobstrct=F
    ##!!
  elif [[ $wamresol = swamp ]] ; then
    # swamp setup
    xdella=0.0500000
    amosop=-2.50000
    amonop=2.50000
    amowep=0.00000
    amoeap=5.0000
    deptha=1.0
    iper=0
    irgg=0
    ##!!
    cldomain=s
    llobstrct=F
    ##!!
  elif [[ $wamresol = standrewbay ]] ; then
    # St Andrew Bay setup
    xdella=0.0010000
    amosop=30.02000
    amonop=30.18000
    amowep=-85.74700
    amoeap=-85.6450
    deptha=0.1
    iper=0
    irgg=0
    ##!!
    cldomain=m
    llobstrct=F
    fr1=0.3
    ##!!
  elif [[ $wamresol = gm20 ]] ; then
    # Gulf of Mexico
    xdella=0.2000000
    amosop=18.00000
    amonop=31.00000
    amowep=-98.00000
    amoeap=-80.0000
    iper=0
    irgg=0
    ##!!
    cldomain=m
    llobstrct=F
    ##!!
  elif [[ ${wamresol} = custom ]] ; then
    xdella=3.0
    irgg=0
    amosop=-78.0
    amonop=81.0
    amowep=0.0
    amoeap=357.0
    iper=1
  else
    echo "$0: ERROR -  wamresol=$wamresol.  Invalid resolution."
    echo "$0: ERROR -  EXIT"
    exit
  fi

  set +eu 

  #======================================================================
  #
  #  End of wave_set_config 
  #
  #======================================================================
  set -e; return 0
}

