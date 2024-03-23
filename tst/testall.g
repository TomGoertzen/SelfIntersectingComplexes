#############################################################################
##  testall.g
#############################################################################
##
##  This file runs package tests.
##
#############################################################################
##
##  This file is part of the SelfIntersectingComplexes package.
##
##  This file's authors include Christian Amend and Tom Goertzen.
##
##  Please refer to the COPYRIGHT file for details.
##
##  SPDX-License-Identifier: GPL-2.0-or-later
##
#############################################################################

LoadPackage("SelfIntersectingComplexes");

TestDirectory(DirectoriesPackageLibrary( "SelfIntersectingComplexes", "tst/Files" ),
  rec(exitGAP := true));




FORCE_QUIT_GAP(1); # if we ever get here, there was an error
