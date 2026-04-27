using BolometricCorrections: MISTChemistry, MISTChemistryv1, MISTBCGrid, MISTBCGridv1, MISTBCTable, MISTBCTablev1
import Logging
using Test

Logging.with_logger(Logging.ConsoleLogger(stderr, Logging.Error)) do
    # inner test to check that the deprecation warning is thrown, that returns the result of the test, 
    # which should be true if MISTChemistry() is an alias for MISTChemistryv1(), so we wrap in a test
    # @test(@test_warn "deprecated" MISTChemistryv1() == MISTChemistry()) # Test that MISTChemistry() is not v2

    # grid = MISTBCGridv1("GALEX")
    # @test(@test_warn "deprecated" grid.table == MISTBCGrid("GALEX").table) # Test that MISTBCGrid() is an alias for MISTBCGridv1()
    # table = grid(-1.0, 0.0)
    # @test(@test_warn "deprecated" table.itp == MISTBCTable(grid, -1, 0.0).itp)

    @test MISTChemistryv1() == MISTChemistry()

    grid = MISTBCGridv1("GALEX")
    @test grid.table == MISTBCGrid("GALEX").table
    table = grid(-1.0, 0.0)
    @test table.itp == MISTBCTable(grid, -1, 0.0).itp
end