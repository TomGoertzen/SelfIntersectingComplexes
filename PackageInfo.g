#############################################################################
##  PackageInfo.g
#############################################################################
##
##  This file contains package meta data. For additional information on
##  the meaning and correct usage of these fields, please consult the
##  manual of the "Example" package as well as the comments in its
##  PackageInfo.g file.
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


SetPackageInfo( rec(

PackageName := "SelfIntersectingComplexes",
Subtitle := "provides algorithms for retriangulating and computing outer hulls of self-intersecting complexes.",
Version := "1.0dev",
Date := "11/03/2024", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Christian",
    LastName := "Amend",
    #WWWHome := TODO,
    Email := "christian.amend@rwth-aachen.de",
    #PostalAddress := TODO,
    #Place := TODO,
    #Institution := TODO,
  ),
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Tom",
    LastName := "Goertzen",
    #WWWHome := TODO,
    Email := "tom.goertzen@rwth-aachen.de",
    #PostalAddress := TODO,
    #Place := TODO,
    #Institution := TODO,
  ),
],



SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/TomGoertzen/SelfIntersectingComplexes",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
#PackageWWWHome := "https://TomGoertzen.github.io/SelfIntersectingComplexes/",
#PackageInfoURL := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
#README_URL     := Concatenation( ~.PackageWWWHome, "README.md" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=  "",

PackageDoc := rec(
  BookName  := "SelfIntersectingComplexes",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "provides algorithms for retriangulating and computing outer hulls of self-intersecting complexes.",
),

Dependencies := rec(
  GAP := ">= 4.12",
  NeededOtherPackages := [["SimplicialSurfaces",">=0.6"]],
  SuggestedOtherPackages := [["GAPic",">=0.1"]],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

#Keywords := [ "TODO" ],

));


