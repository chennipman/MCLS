MCLS
====

TODO

Build
=====

To build the main executable, run

    make

or

    make MCLS

If you want to build in parallel, say on 4 threads, run

    make -j4

You can run all unit tests by calling

    make test

or just one, e.g.

    make test_momentum_predictor

All unit tests in the `unittest` directory are automatically available with the
`test_` prefix.
