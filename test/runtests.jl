using GlyphEctoplasm
using Test
using Mustache
using Colors
using FixedPointNumbers: N0f8
using Random
using GlyphEctoplasm.PNGFiles

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

    @testset "PNG quantization" begin
        Random.seed!(1)
        dir = mktempdir()
        try
            # A noisy truecolor image (many distinct colors) is large as RGB and
            # shrinks substantially once reduced to a small palette.
            img = rand(RGB{N0f8}, 128, 128)
            p = joinpath(dir, "a.png")
            PNGFiles.save(p, img)
            before = filesize(p)

            ncolors = quantize_png!(p; ncolors=16)
            @test 0 < ncolors <= 16            # palette built, within budget
            @test filesize(p) < before          # actually smaller

            # Output is still a valid PNG of identical dimensions.
            out = PNGFiles.load(p)
            @test size(out) == (128, 128)

            # Re-running is safe (idempotent) and stays within budget.
            @test 0 < quantize_png!(p; ncolors=16) <= 16

            # Directory sweep reports a real reduction.
            PNGFiles.save(joinpath(dir, "b.png"), rand(RGB{N0f8}, 64, 64))
            r = optimize_pngs!(dir; verbose=false)
            @test r.files == 2
            @test r.bytes_after <= r.bytes_before

            # Graceful failure: non-PNG input and missing file return -1, no throw.
            txt = joinpath(dir, "not_an_image.png")
            write(txt, "this is not a png")
            @test quantize_png!(txt) == -1
            @test quantize_png!(joinpath(dir, "missing.png")) == -1
        finally
            rm(dir; recursive=true, force=true)
        end
    end

    @testset "bin_value" begin
        bv = GlyphEctoplasm.bin_value
        @test bv(0.0, 0.0, 10.0, 10) == 0
        @test bv(10.0, 0.0, 10.0, 10) == 9      # top value clamps into last bin
        @test bv(5.0, 0.0, 10.0, 10) == 5
        @test bv(2.5, 0.0, 10.0, 10) == 2
        @test bv(NaN, 0.0, 10.0, 10) == -1      # missing -> sentinel (sorts last)
        @test bv(3.0, 2.0, 2.0, 10) == 0        # degenerate range -> bin 0
    end

    @testset "binned sort ordering" begin
        # (label, group_order, shapley_median, cluster_median, count)
        # group_order: 0 = single_region, 2 = 2_regions, 3 = 3_regions, ...
        motifs = [
            ("single",        0, 0.1,  5.0, 10),
            ("r3_first",      3, 0.9,  1.0, 50),   # higher group sorts after lower group
            ("r2_hiShapHiC",  2, 0.95, 9.0,  5),
            ("r2_hiShapLoC",  2, 0.96, 0.5,  7),
            ("r2_hiShapLoC2", 2, 0.97, 0.4, 99),   # same shap+cluster bin, higher count
            ("r2_loShap",     2, -0.9, 8.0, 100),
            ("r2_noCluster",  2, 0.5,  NaN,  3),
        ]
        bc = 10
        slo, shi = extrema(m[3] for m in motifs if !isnan(m[3]))
        clo, chi = extrema(m[4] for m in motifs if !isnan(m[4]))
        key(m) = (m[2],                                              # group order primary
                  -GlyphEctoplasm.bin_value(m[3], slo, shi, bc),
                  -GlyphEctoplasm.bin_value(m[4], clo, chi, bc),
                  -m[5])
        labels = [m[1] for m in sort(motifs, by=key)]

        @test labels[1] == "single"             # single_region group leads
        # all 2_regions come before the 3_regions group, regardless of influence
        @test findfirst(==("r3_first"), labels) > findfirst(==("r2_loShap"), labels)
        # within 2_regions: high shap bin first, then high cluster bin
        @test labels[2] == "r2_hiShapHiC"       # higher cluster bin first within shap bin
        @test labels[3] == "r2_hiShapLoC2"      # count breaks tie (99 > 7) within same bins
        @test labels[4] == "r2_hiShapLoC"
        @test labels[end] == "r3_first"         # last group (3_regions) sorts last
    end

    @testset "cluster median display text" begin
        paths = GlyphEctoplasm.build_motif_paths("x", mktempdir(), "singletons")
        # conv-style result (no cluster_median field) must not crash or show it
        t1 = build_metadata_texts(nothing, paths, 0.5, 10; show_meme_and_csv=false,
                nnd_result=(; k=5, obs_mNND=0.12, p_value=0.03))
        @test !occursin("cluster median", t1[end])
        # mut-style result shows the value
        t2 = build_metadata_texts(nothing, paths, 0.5, 10; show_meme_and_csv=false,
                nnd_result=(; k=5, obs_mNND=0.12, p_value=0.03, cluster_median=1.234))
        @test occursin("cluster median: <strong>1.234", t2[end])
        # NaN is omitted
        t3 = build_metadata_texts(nothing, paths, 0.5, 10; show_meme_and_csv=false,
                nnd_result=(; k=5, obs_mNND=0.12, p_value=0.03, cluster_median=NaN))
        @test !occursin("cluster median", t3[end])
    end

    # ── location statistic (companion to the NND tightness test) ────────────

    @testset "location_z_1d" begin
        lz = GlyphEctoplasm.location_z_1d

        # Hand-checkable closed form: background 1..5 has mean 3 and sample sd
        # sqrt(2.5); carriers {1,2,3} have mean 2 and n=3.
        bg_small = [1.0, 2.0, 3.0, 4.0, 5.0]
        @test lz([1, 2, 3], bg_small) ≈ (2.0 - 3.0) / (sqrt(2.5) / sqrt(3))

        Random.seed!(20250820)
        bg = randn(4000)

        # Strongly shifted carriers -> large |z|, sign follows the shift.
        up = collect(1:150)
        bg_up = copy(bg); bg_up[up] .+= 3.0
        z_up = lz(up, bg_up)
        @test z_up > 10
        bg_dn = copy(bg); bg_dn[up] .-= 3.0
        z_dn = lz(up, bg_dn)
        @test z_dn < -10

        # A random subset of the same size is a draw from the null -> |z| ~ N(0,1).
        rand_pos = Random.randperm(length(bg))[1:150]
        @test abs(lz(rand_pos, bg)) < 3

        # Degenerate cases return the NaN sentinel instead of throwing.
        @test isnan(lz(Int[], bg))                    # no carriers
        @test isnan(lz([7], bg))                      # one carrier
        @test isnan(lz([1, 2, 3], fill(2.0, 20)))     # zero-variance background
        @test isnan(lz([1, 2, 3], nothing))           # missing labels
        @test isnan(lz(nothing, bg))                  # missing carrier list
        @test isnan(lz([1, 2], [4.0]))                # background too small

        # `missing` / non-finite labels are skipped, not propagated.
        bg_m = Vector{Union{Missing,Float64}}(bg_small)
        bg_m[5] = missing
        @test lz([1, 2, 3], bg_m) ≈ (2.0 - 2.5) / (sqrt(5/3) / sqrt(3))
        @test isnan(lz([1, 5], bg_m))                 # only 1 usable carrier left
        bg_nan = [1.0, 2.0, 3.0, 4.0, NaN, Inf]
        @test lz([1, 2, 3], bg_nan) ≈ lz([1, 2, 3], [1.0, 2.0, 3.0, 4.0])

        # Positional indices outside the background are ignored.
        @test lz([1, 2, 3, 999], bg_small) ≈ lz([1, 2, 3], bg_small)
    end

    @testset "location z display text" begin
        paths = GlyphEctoplasm.build_motif_paths("x", mktempdir(), "singletons")
        res = (; k=5, obs_mNND=0.12, p_value=0.03, cluster_median=1.234, location_z=-6.7)

        # Legacy call (no `report_location_z`): unchanged, even though the
        # nnd_result now carries the field.
        legacy = build_metadata_texts(nothing, paths, 0.5, 10; show_meme_and_csv=false,
                nnd_result=res)
        @test !occursin("location z", legacy[end])
        @test legacy == build_metadata_texts(nothing, paths, 0.5, 10; show_meme_and_csv=false,
                nnd_result=(; k=5, obs_mNND=0.12, p_value=0.03, cluster_median=1.234))

        # Opted in: shown beside the NND p-value, sign preserved.
        shown = build_metadata_texts(nothing, paths, 0.5, 10; show_meme_and_csv=false,
                nnd_result=res, report_location_z=true)
        @test occursin("p-value (NND)", shown[end])
        @test occursin("location z: <strong>-6.70", shown[end])

        # NaN sentinel is omitted rather than printed.
        nan_z = build_metadata_texts(nothing, paths, 0.5, 10; show_meme_and_csv=false,
                nnd_result=(; k=5, obs_mNND=0.12, p_value=0.03, cluster_median=1.234,
                              location_z=NaN), report_location_z=true)
        @test !occursin("location z", nan_z[end])

        # A conv-style result predating the field: opting in is a no-op, no throw.
        old_shape = build_metadata_texts(nothing, paths, 0.5, 10; show_meme_and_csv=false,
                nnd_result=(; k=5, obs_mNND=0.12, p_value=0.03), report_location_z=true)
        @test !occursin("location z", old_shape[end])
    end

    @testset "TopMoverEntry location_z back-compat" begin
        args20 = (0.5, 1.0, 0.2, 30, "motif 1", "Singleton Motifs", "", 0.01,
                  true, "", GlyphEctoplasm.WildTypeRegion[], "a.png", "", "",
                  "i.png", "k.png", ["t"], String[], String[], Vector{String}[])

        # 20-argument (pre-change) construction still works, marks the stat absent.
        legacy = TopMoverEntry(args20...)
        @test isnan(legacy.location_z)
        # 21-argument construction carries the value; every other field matches.
        withz = TopMoverEntry(args20..., -8.25)
        @test withz.location_z == -8.25
        @test all(getfield(legacy, f) == getfield(withz, f)
                  for f in fieldnames(TopMoverEntry) if f !== :location_z)

        # CSV schema: absent by default, so existing readers see the same table.
        base = GlyphEctoplasm.top_movers_dataframe([withz], TopMoverEntry[])
        @test !("location_z" in names(base))
        @test names(base)[12:13] == ["cluster_nnd", "nnd_pvalue"]

        # Opted in: one extra column, immediately after `nnd_pvalue`.
        withcol = GlyphEctoplasm.top_movers_dataframe([withz], TopMoverEntry[];
                                                      report_location_z=true)
        @test names(withcol) == vcat(names(base)[1:13], "location_z", names(base)[14:end])
        @test withcol.location_z == [-8.25]
        @test isnan(GlyphEctoplasm.top_movers_dataframe([legacy], TopMoverEntry[];
                        report_location_z=true).location_z[1])
    end

end
