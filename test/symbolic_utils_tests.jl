using QuasiGeostrophy
using Test
using SymbolicUtils
import SymbolicUtils: Chain, Postwalk
import SymbolicUtils: Sym, Term, istree, operation, arguments, to_symbolic, Fixpoint

include(pwd() * "/test/test_utils.jl")

@wrapper a=1 b=2

@testset "symbolic utils swapping" begin
    println("for $a = ", a.data, ", $b = ", b.data)
    println("changing a+b to a*b")
    c = a+b
    check_c = a*b
    r1 = @acrule ~ra + ~rb => ~ra * ~rb
    new_c = Fixpoint(Postwalk(Chain([r1])))(c)
    @test compute(check_c) == compute(new_c)
    println("changing sin(a+b) * sin(a+a+b) to sin(a*b) * sin(a*a*b)")
    c = sin(a+b) * sin(a+a+b)
    check_c = sin(a*b) * sin(a*a*b)
    new_c = Fixpoint(Postwalk(Chain([r1])))(c)
    @test compute(check_c) == compute(new_c)
    println("Checking a*a = a^2 with rewrite rules")
    r2 = @acrule ~ra * ~ra => Exponentiation(~ra , 2)
    c = a * a
    new_c = Fixpoint(Postwalk(Chain([r2])))(c)
    @test (new_c != c) * (compute(c) == compute(new_c))
end