#! /usr/bin/env python

from __future__ import print_function

import re
import os


def parse( filename, funcdefs ):

    with open( filename, 'r' ) as file:
        file = file.read()
    # remove all comments
    file, n = re.subn( r'//.*$', lambda m: '', file, flags = re.MULTILINE )
    file, n = re.subn( r'/[*].*?[*]/', lambda m: '', file, flags = re.DOTALL )
    # find EXPORT and write to `funcdefs`
    for match in re.findall( r'EXPORT .*?\(.*?\)', file, re.DOTALL ):
        match = match[7:] + ';'
        match, n = re.subn( r'\s+', lambda m: ' ', match )
        match, n = re.subn( r'\s*,\s?', lambda m: ', ', match )
        match, n = re.subn( r'\s*\)\s*', lambda m: ' )', match )
        match, n = re.subn( r'\s*\(\s*', lambda m: '( ', match )
        match, n = re.subn( r'([(,] ?)(string|ofstream)', lambda m: m.group(1) + 'std::' + m.group(2), match ) # )
        print( match, file = funcdefs )


if __name__ == '__main__':

    with open( '../objects/funcdefs.h', 'w' ) as funcdefs:
        print( '#pragma once', file = funcdefs )
        print( '#define EXPORT', file = funcdefs )
        print( '#include "../src/headers/header_clsses.h"', file = funcdefs )
        print( '#include "../src/headers/array.h"', file = funcdefs )
        for root, dirs, files in os.walk( '.' ):
            for file in files:
                if not file.endswith( '.cpp' ):
                    continue
                parse( os.path.join( root, file ), funcdefs )
