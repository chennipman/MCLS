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

    import sys

    source_files = []
    output = None

    args = iter( sys.argv )
    next( args ) # skip program name
    for arg in args:
        if arg == '--':
            break
        elif arg == '--output':
            output = next( args )
        else:
            source_files.append( arg )
    source_files.extend( args )

    if output is None:
        print( 'ERROR: `--output OUTPUT` argument missing' )
        raise SystemExit( 1 )

    with open( output, 'w' ) as funcdefs:
        print( '#pragma once', file = funcdefs )
        print( '#define EXPORT', file = funcdefs )
        path_to_headers = os.path.relpath( 'src/headers', os.path.dirname( output ) )
        print( '#include "{}/header_clsses.h"'.format( path_to_headers ), file = funcdefs )
        print( '#include "{}/array.h"'.format( path_to_headers ), file = funcdefs )
        for file in source_files:
            if not file.endswith( '.cpp' ):
                continue
            parse( file, funcdefs )
