#!/usr/bin/python
# -*- coding: UTF-8 -*-
from __future__ import absolute_import

from re import match,findall,compile
from math import log,e,inf
from decimal import Decimal as dec
from scipy import stats
from progressbar import ProgressBar,ETA,Bar,Percentage
from sys import exit
import copy
from .candris import open_file,count_mutation,prepare_fas_cds_file,cncs,candris,H_test,two_component_cncs
from pkg_resources import resource_filename as data_path
__all__ = ["candris"]