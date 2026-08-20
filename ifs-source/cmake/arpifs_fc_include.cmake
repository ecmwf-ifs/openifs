# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

# Direct arpifs source add-back list removed for OpenIFS minimisation.
# Keep pattern-based OpenIFS support source additions below.

ecbuild_list_add_pattern(LIST arpifs.${PREC}_src GLOB
# need to add the local copy of emos libs and forecast-only dummies 
# so that arpifs can build with both forecast-only and openifs-only
  openifs/emos/*
  openifs/dummy_fc_only/*
  openifs/dummy_ifsobs/*
  openifs/var/*
)
