using BolometricCorrections: MISTChemistry, MISTv1Chemistry, MISTBCGrid, MISTv1BCGrid, MISTBCTable, MISTv1BCTable
import Logging
using Test

Logging.with_logger(Logging.ConsoleLogger(stderr, Logging.Error)) do
    # inner test to check that the deprecation warning is thrown, that returns the result of the test, 
    # which should be true if MISTChemistry() is an alias for MISTv1Chemistry(), so we wrap in a test
    # @test(@test_warn "deprecated" MISTv1Chemistry() == MISTChemistry()) # Test that MISTChemistry() is not v2

    # grid = MISTv1BCGrid("GALEX")
    # @test(@test_warn "deprecated" grid.table == MISTBCGrid("GALEX").table) # Test that MISTBCGrid() is an alias for MISTv1BCGrid()
    # table = grid(-1.0, 0.0)
    # @test(@test_warn "deprecated" table.itp == MISTBCTable(grid, -1, 0.0).itp)

    @test MISTv1Chemistry() == MISTChemistry()

    grid = MISTv1BCGrid("GALEX")
    @test grid.table == MISTBCGrid("GALEX").table
    table = grid(-1.0, 0.0)
    @test table.itp == MISTBCTable(grid, -1, 0.0).itp
end