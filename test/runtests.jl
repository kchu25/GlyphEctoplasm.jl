using GlyphEctoplasm
using Test
using Mustache

@testset "GlyphEctoplasm.jl" begin
    
    @testset "init_json_dict" begin
        # Empty initialization
        d = init_json_dict()
        @test d isa Dict{String, Dict{String, GlyphEctoplasm.JSONValue}}
        @test isempty(d)
        
        # With pre-populated keys
        d2 = init_json_dict(keys=["singletons", "pairs"])
        @test haskey(d2, "singletons")
        @test haskey(d2, "pairs")
        @test isempty(d2["singletons"])
    end
    
    @testset "init_dict_for_html_render" begin
        html_dict = init_dict_for_html_render()
        
        # Check expected keys exist
        @test haskey(html_dict, "div_img_id")
        @test haskey(html_dict, "i")
        @test haskey(html_dict, "img_src")
        @test haskey(html_dict, "p_id1_default")
        @test haskey(html_dict, "p_id7_default")
        @test haskey(html_dict, "group_id")
        @test haskey(html_dict, "button_text")
        @test haskey(html_dict, "has_interaction")
        
        # All values should be empty vectors
        @test all(isempty(v) for v in values(html_dict))
        @test all(v isa Vector{String} for v in values(html_dict))
    end
    
    @testset "ConvMotifConfig construction" begin
        # Test struct definition exists
        @test ConvMotifConfig <: Any
        @test fieldnames(ConvMotifConfig) == (
            :data, :filter_len, :float_type, :dpi, :alpha, :use_rna, :xlim, :save_path
        )
    end
    
    @testset "MutationRegionConfig construction" begin
        # Test struct definition
        @test GlyphEctoplasm.MutationRegionConfig <: Any
        @test :dpi in fieldnames(GlyphEctoplasm.MutationRegionConfig)
        @test :use_rna in fieldnames(GlyphEctoplasm.MutationRegionConfig)
        @test :save_path in fieldnames(GlyphEctoplasm.MutationRegionConfig)
    end
    
    @testset "HTML templates exist" begin
        # Templates are Mustache.MustacheTokens
        @test html_template_unified isa Mustache.MustacheTokens
        @test html_template_generalization isa Mustache.MustacheTokens
        @test script_template isa Mustache.MustacheTokens
        @test template_css isa Mustache.MustacheTokens
    end
    
    @testset "Constants defined" begin
        # Check nucleotide mappings
        @test GlyphEctoplasm._ind2dna_str_[1] == 'A'
        @test GlyphEctoplasm._ind2dna_str_[4] == 'T'
        @test GlyphEctoplasm._ind2dna_str_rna[4] == 'U'
        
        # Check placeholder
        @test GlyphEctoplasm._placeholder_char_ == '-'
    end

end
