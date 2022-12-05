##########################################################
# Copyright (c) 2022 datavisyn GmbH, http://datavisyn.io
#
# This file is property of datavisyn.
# Code and any other files associated with this project
# may not be copied and/or distributed without permission.
#
# Proprietary and confidential. No warranty.
#
##########################################################
from tdp_core.dbview import DBConnector


def create():
    return DBConnector(views=dict())
