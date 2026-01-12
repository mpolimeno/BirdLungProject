@testset "COMSOL data parsing" begin

    # Test parameters
    t_initial = 0.0
    time_step = 0.02
    t_final   = 5.86
    nNodes    = 5
    test_time_length = size(collect(t_initial:time_step:t_final),1)

    @testset "x-coordinate parsing" begin

        # Values of x_coordinate (first 5 points) at t=0, t=0.02 and t=5.86 taken from line C
        x_t0 = [
            0.1875,
            0.1875,
            0.1875,
            0.1875,
            0.1875
        ]

        x_t02 = [
            0.1875,
            0.1875,
            0.1875,
            0.1875,
            0.1875
        ]

        x_t5p86 = [
            0.1875,
            0.1875,
            0.1875,
            0.1875,
            0.1875
        ]

        # create a temporary file with the listed values for x
        mktempdir() do dir
            filename = joinpath(dir, "comsol_test_x_coordinate.txt")

            open(filename, "w") do io
                write(io, """
                % Data
                % x (m) @ t=0
                $(join(x_t0, " "))

                % Data
                % x (m) @ t=0.02
                $(join(x_t02, " "))

                % Data
                % x (m) @ t=5.86
                $(join(x_t5p86, " "))
                """)
            end

            # call function
            x_coordinate = DataParsing.extract_xcoordinatedata_from_COMSOL(
                filename,
                t_initial,
                time_step,
                t_final,
                nNodes
            )

            # check if listed values and the one from the files agree
            @test size(x_coordinate) == (nNodes, test_time_length)
            @test isapprox(x_coordinate[:,1],x_t0)
            @test isapprox(x_coordinate[:,2],x_t02)
            @test isapprox(x_coordinate[:,test_time_length],x_t5p86)  
        end
    end

    @testset "y_coordinate parsing" begin
        # Values of y_coordinate (first 5 points) at t=0, t=0.02 and t=5.86 taken from line C
        y_t0 = [
            -0.00497499999999999,
            -0.00492499999999999,
            -0.00487499999999999,
            -0.004824999999999991,
            -0.00477499999999999
        ]

        y_t02 = [
            -0.00497499999999999,
            -0.00492499999999999,
            -0.00487499999999999,     
            -0.004824999999999991,
            -0.00477499999999999
        ]

        y_t5p86 = [
            -0.00497499999999999,
            -0.00492499999999999,
            -0.00487499999999999,
            -0.004824999999999991,
            -0.00477499999999999
        ]

        # create a temporary file with the listed values for y
        mktempdir() do dir
            filename = joinpath(dir, "comsol_test_y_coordinate.txt")

            open(filename, "w") do io
                write(io, """
                % Data
                % y (m) @ t=0
                $(join(y_t0, " "))

                % Data
                % y (m) @ t=0.02
                $(join(y_t02, " "))

                % Data
                % y (m) @ t=5.86
                $(join(y_t5p86, " "))
                """)
            end

            # call function
            y_coordinate = DataParsing.extract_ycoordinatedata_from_COMSOL(
                filename,
                t_initial,
                time_step,
                t_final,
                nNodes
            )

            # check if listed values and the one from the files agree
            @test size(y_coordinate) == (nNodes, test_time_length)
            @test isapprox(y_coordinate[:,1],y_t0)
            @test isapprox(y_coordinate[:,2],y_t02)
            @test isapprox(y_coordinate[:,test_time_length],y_t5p86)  
        end
    end

    @testset "horizontal velocity parsing" begin

        # Values of u (first 5 points) at t=0, t=0.02 and t=5.86 taken from line C
        u_t0 = [
            -1.244093008126981E-7,
            -3.5526567789660583E-7,
            -5.62172422258546E-7,
            -7.451295338985086E-7,
            -9.041370128165054E-7,
        ]

        u_t02 = [
            -3.553451898117511E-5,
            -1.0197458909379784E-4,
            -1.622427020734555E-4,
            -2.1633885792014514E-4,
            -2.6426305663387027E-4
        ]

        u_t5p86 = [
            -1.9974638636147816E-5,
            -4.881703935446194E-5,
            -6.285027133413621E-5,
            -6.207433457517028E-5,
            -4.648922907756411E-5
        ]

        # create a temporary file with the listed values for u
        mktempdir() do dir
            filename = joinpath(dir, "comsol_test_u_velocity.txt")

            open(filename, "w") do io
                write(io, """
                % Data
                % u (m/s) @ t=0
                $(join(u_t0, " "))

                % Data
                % u (m/s) @ t=0.02
                $(join(u_t02, " "))

                % Data
                % u (m/s) @ t=5.86
                $(join(u_t5p86, " "))
                """)
            end

            # call function
            u_velocity = DataParsing.extract_horizontalvelocity_data_from_COMSOL(
                filename,
                t_initial,
                time_step,
                t_final,
                nNodes
            )

            # check if listed values and the one from the files agree
            @test size(u_velocity) == (nNodes, test_time_length)
            @test isapprox(u_velocity[:,1],u_t0)
            @test isapprox(u_velocity[:,2],u_t02)
            @test isapprox(u_velocity[:,test_time_length],u_t5p86)  
        end
    end
    
    @testset "vertical velocity parsing" begin

        # Values of u (first 5 points) at t=0, t=0.02 and t=5.86 taken from line C
        v_t0 = [
            -1.6216119084450657E-10,
            -4.72450835053515E-10,
            -7.640301626225398E-10,
            -1.0368991735515664E-9,
            -1.2910578678406129E-9
        ]

        v_t02 = [
            -9.37875247593456E-8,
            -2.8136606460193274E-7,
            -4.6894925820972756E-7,
            -6.565371055827206E-7,
            -8.441296067209242E-7
        ]

        v_t5p86 = [
            -1.026838877180644E-7,
            -3.094535967071683E-7,
            -5.180925504335868E-7,
            -7.286007488973091E-7,
            -9.409781920983491E-7
        ]

        # create a temporary file with the listed values for v
        mktempdir() do dir
            filename = joinpath(dir, "comsol_test_v_velocity.txt")

            open(filename, "w") do io
                write(io, """
                % Data
                % v (m/s) @ t=0
                $(join(v_t0, " "))

                % Data
                % v (m/s) @ t=0.02
                $(join(v_t02, " "))

                % Data
                % v (m/s) @ t=5.86
                $(join(v_t5p86, " "))
                """)
            end

            # call function
            v_velocity = DataParsing.extract_verticalvelocity_data_from_COMSOL(
                filename,
                t_initial,
                time_step,
                t_final,
                nNodes
            )

            # check if listed values and the one from the files agree
            @test size(v_velocity) == (nNodes, test_time_length)
            @test isapprox(v_velocity[:,1],v_t0)
            @test isapprox(v_velocity[:,2],v_t02)
            @test isapprox(v_velocity[:,test_time_length],v_t5p86)  
        end
    end

    @testset "pressure parsing" begin

        # Values of pressure (first 5 points) at t=0, t=0.02 and t=5.86 taken from line C
        p_t0 = [
            7.261650049132611,
            7.261650387837785,
            7.26165072654296,
            7.261651065248136,
            7.261651403953311
        ]

        p_t02 = [
            18.70010308337058,
            18.700103014188983,
            18.70010294500739,
            18.7001028758258,
            18.700102806644203
        ]

        p_t5p86 = [
            16.253604862340325,
            16.253604899969716,
            16.253604937599103,
            16.253604975228495,
            16.253605012857882
        ]

        # create a temporary file with the listed values for the pressure
        mktempdir() do dir
            filename = joinpath(dir, "comsol_test_pressure.txt")

            open(filename, "w") do io
                write(io, """
                % Data
                % p (Pa) @ t=0
                $(join(p_t0, " "))

                % Data
                % p (Pa) @ t=0.02
                $(join(p_t02, " "))

                % Data
                % p (Pa) @ t=5.86
                $(join(p_t5p86, " "))
                """)
            end

            # call function
            pressure = DataParsing.extract_pressuredata_from_COMSOL(
                filename,
                t_initial,
                time_step,
                t_final,
                nNodes
            )

            # check if listed values and the one from the files agree
            @test size(pressure) == (nNodes, test_time_length)
            @test isapprox(pressure[:,1],p_t0)
            @test isapprox(pressure[:,2],p_t02)
            @test isapprox(pressure[:,test_time_length],p_t5p86)  
        end
    end

    @testset "horizontal acceleration parsing" begin

        # Values of ut (first 5 points) at t=0, t=0.02 and t=5.86 taken from line C
        ut_t0 = [
            -7.976131648791424E-4,
            -0.00227771539403278,
            -0.0036043188223803314,
            -0.004777423449921731,   
            -0.005797029276657055
        ]

        ut_t02 = [
            -0.0016627385307387737,    
            -0.004785239992899596,
            -0.007637107322638351,
            -0.010218340519954896,   
            -0.012528939584849403
        ]

        ut_t5p86 = [
            -8.786799378591801E-4,
            -0.002571798885450865,    
            -0.004179263262207108,    
            -0.00570107306812783,     
            -0.007137228303213129
        ]

        # create a temporary file with the listed values for ut
        mktempdir() do dir
            filename = joinpath(dir, "comsol_test_u_acceleration.txt")

            open(filename, "w") do io
                write(io, """
                % Data
                % ut (m/s^2) @ t=0
                $(join(ut_t0, " "))

                % Data
                % ut (m/s^2) @ t=0.02
                $(join(ut_t02, " "))

                % Data
                % ut (m/s^2) @ t=5.86
                $(join(ut_t5p86, " "))
                """)
            end

            # call function
            u_acceleration = DataParsing.extract_horizontalacceleration_data_from_COMSOL(
                filename,
                t_initial,
                time_step,
                t_final,
                nNodes
            )

            # check if listed values and the one from the files agree
            @test size(u_acceleration) == (nNodes, test_time_length)
            @test isapprox(u_acceleration[:,1],ut_t0)
            @test isapprox(u_acceleration[:,2],ut_t02)
            @test isapprox(u_acceleration[:,test_time_length],ut_t5p86)  
        end
    end

    @testset "vertical acceleration parsing" begin

        # Values of vt (first 5 points) at t=0, t=0.02 and t=5.86 taken from line C
        vt_t0 = [
            -1.0656009152716386E-6,    
            -3.1082254648596057E-6,   
            -5.032746973173974E-6,    
            -6.839165440214649E-6,
            -8.527480865981746E-6
        ]

        vt_t02 = [
            -3.1492805267828837E-6,
            -9.457935710550285E-6,
            -1.5780049734586967E-5,   
            -2.2115622598892618E-5,   
            -2.846465430346764E-5    
        ]

        vt_t5p86 = [
            -1.2948530721136508E-8,
            -4.025961059992342E-8,
            -6.9456048394064E-8,
            -1.0053784410355675E-7,
            -1.335049977284037E-7
        ]

        # create a temporary file with the listed values for vt
        mktempdir() do dir
            filename = joinpath(dir, "comsol_test_v_acceleration.txt")

            open(filename, "w") do io
                write(io, """
                % Data
                % vt (m/s^2) @ t=0
                $(join(vt_t0, " "))

                % Data
                % vt (m/s^2) @ t=0.02
                $(join(vt_t02, " "))

                % Data
                % vt (m/s^2) @ t=5.86
                $(join(vt_t5p86, " "))
                """)
            end

            # call function
            v_acceleration = DataParsing.extract_verticalacceleration_data_from_COMSOL(
                filename,
                t_initial,
                time_step,
                t_final,
                nNodes
            )

            # check if listed values and the one from the files agree
            @test size(v_acceleration) == (nNodes, test_time_length)
            @test isapprox(v_acceleration[:,1],vt_t0)
            @test isapprox(v_acceleration[:,2],vt_t02)
            @test isapprox(v_acceleration[:,test_time_length],vt_t5p86)  
        end
    end

end
