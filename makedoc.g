#############################################################################
##  makedoc.g
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


if fail = LoadPackage("AutoDoc", "2018.02.14") then
    Error("AutoDoc version 2018.02.14 or newer is required.");
fi;

AutoDoc( rec( scaffold := rec(
        includes := [
            "intro.xml",
            "functions.xml",
            "license.xml",
            ],
        ),
        extract_examples := true,
        autodoc := true ) );

Exec("dev/processTests.sh");
