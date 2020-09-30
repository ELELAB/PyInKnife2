#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    dask_patches.py
#
#    Dask monkey patches.
#
#    Copyright (C) 2020 Valentina Sora 
#                       <sora.valentina1@gmail.com>
#                       Juan Salamanca Viloria 
#                       <juan.salamanca.viloria@gmail.com> 
#                       Elena Papaleo
#                       <elenap@cancer.dk>
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program. 
#    If not, see <http://www.gnu.org/licenses/>.

import dask
import distributed.utils
import logging


# to address a bug that resets the distributed.worker
# logger to WARNING level when a task is launched on
# the worker, no matter what the configuration was
def reset_worker_logger():
    """Utility function to reset a Dask logger handlers
    and level to desired values.
    """

    # new level
    NEWLEVEL = logging.INFO
    # get the logger
    logger = logging.getLogger("distributed.worker")
    # define the handlers to keep
    htokeep = [h for h in logger.handlers if type(h).__name__ == \
               distributed.utils.DequeHandler.__name__]
    # remove all the handlers
    for h in logger.handlers:
        logger.removeHandler(h)
    # add the handlers to keep
    for h in htokeep:
        # set the new level
        h.setLevel(NEWLEVEL)
        # add the handler to the logger
        logger.addHandler(h)
    # reset the logger level to the new level
    logger.setLevel(NEWLEVEL)
    # return the new logger
    return logger