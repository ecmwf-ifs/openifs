# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction
# 
# (C) Copyright 1989- Meteo-France.
# 

"""
A module providing configuration utilities for complex field and
variable data type hierachies.

This module is intended to be used with the Fypp preprocesser to
concisely define ``Variable`` objects from which templated Fortran
type hierarchies can be created through templating.
"""
import yaml

__all__ = ['Variable', 'VariableGroup', 'VariableConfiguration']


class Variable(object):
    """
    Object representing a single variable containing one or more fields.
    """

    def __init__(self, **kwargs):
        self.name = kwargs.get('name')
        self.dim = kwargs.get('dim', 3)
        self.comment = kwargs.get('comment', '')

        # Indicates multi-dimensional array of Variable objects.
        # Note that ``self.array`` normalizes to the variable array
        # rank so scalar variables have ``array==0``.
        self.array = int(kwargs.get('array', False))

    def __repr__(self):
        return self.name


class VariableGroup(object):
    """
    Group of ``Variable`` objects with common dimension and metadata.
    """

    def __init__(self, **kwargs):
        if not kwargs['type'].lower() == 'group':
            raise RuntimeError('Error creating VariableGroup: wrong schema %s' % kwargs)

        # Need a clever way to auto-set attrs
        self.name = kwargs.get('name')
        self.short = kwargs.get('short', self.name)
        self.comment = kwargs.get('comment', '')
        self.dimension = kwargs.get('dimension')
        self.variables = [Variable(dim=self.dimension, **v) for v in kwargs['variables']]

    def __repr__(self):
        return 'Group<%sD>::%s (%s):: %s' % (self.dimension, self.name, self.short, self.variables)


class VariableConfiguration(object):
    """
    Utility class to read a configuration of field variables from .yaml file.
    """

    def __init__(self, filename):
        print("Reading field types configuration from file: %s" % filename)

        with open(filename, 'r') as stream:
            self.schema = yaml.safe_load(stream)

    @property
    def groups(self):
        """
        Get a ``VariableGroup`` from the locally stored variable configuration.
        """
        groups = [VariableGroup(**group) for group in self.schema]
        return {g.name: g for g in groups}
