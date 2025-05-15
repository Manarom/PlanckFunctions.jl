# this testing data was obtained using blackbody calculator 
# https://www.opticsthewebsite.com/OpticsCalculators
module TestingData
export BenchmarkData,benchmark_data,read_temperature_data,temperatures
const citation1 = raw"Berger, Chris G. 'Blackbody Calculator', Optics: The Website, https://www.opticsthewebsite.com/OpticsCalculators . Accessed on 5/14/2025"
"""
This module contains some  `hardcoded` data obtained from $(citation1) 

"""
TestingData

using DelimitedFiles
const spectral_data_folder = joinpath(@__DIR__(),"spectra")

"""   

temperature  - BB temperature, [K]
lower - lower wavelength limit of the band, [μm]
upper - upper wavelength limit of the band, [μm] 
total_radiance - total radiance in the whole spectral range, [W/(m²⋅sr)]
peak_wavelength  - wavelength of Planck function maximum, [μm]
peak_spectral_radiance - spectral radiance at the maximum wavelength, [W/(m²⋅sr⋅μm)]
band_radiance - radiance within the spectral band, [W/(m²⋅sr)]

"""
@kwdef struct BenchmarkData
    temperature # Kelvins
    lower # lower wavelength limit μm
    upper # upper wavelength limit μm
    total_radiance # W/(m²⋅sr)  
    peak_wavelength # μm
    peak_spectral_radiance # W/m²/sr/μm
    band_radiance # W/m²/sr total radiance within the band lower-
    spectrum # 
end
"""
    read_temperature_data(tag)

    Functions to read benchmark data spectrum from the default testing data folder
"""
function read_temperature_data(tag;folder = spectral_data_folder)
    data = readdlm(joinpath(folder,tag),',')
    w = @view data[:,2]
    w.*=1e4 # as far as all data is saved in [W/(cm²⋅sr⋅μm)] we should multiply it
    # by 1e4 to convert units to [W/(m²⋅sr⋅μm)] (deault units in PlanckFunctions)
    return data
end
"""
The data for testing was obtained from

$(citation1)
   

"""
const benchmark_data = Dict(
"300"=>BenchmarkData(temperature = 300.0,
                lower=0.5,upper = 30.0,
                total_radiance=1.46199835e-2*1e4,
                peak_wavelength=9.65923985e+0,
                peak_spectral_radiance =9.95248946e-4*1e4,
                band_radiance=1.30116301e-2*1e4,
                spectrum=read_temperature_data("300.csv") ),
"1000"=>BenchmarkData(temperature = 1000.0,
                lower=0.5,upper = 30.0,
                total_radiance=1.80493624*1e4,
                peak_wavelength=2.89777196,
                peak_spectral_radiance =4.09567468e-1*1e4,
                band_radiance=1.79643710e+0*1e4,
                spectrum=read_temperature_data("1000.csv") ),
"1500"=>BenchmarkData(temperature = 1500.0,
                lower=0.5,upper = 30.0,
                total_radiance=9.13748969e+0*1e4,
                peak_wavelength=1.93184797,
                peak_spectral_radiance =3.11015296e+0*1e4,
                band_radiance=9.12386534*1e4,
                spectrum=read_temperature_data("1500.csv") ),
"2000"=>BenchmarkData(temperature = 2000.0,
                lower=0.5,upper = 30.0,
                total_radiance=	2.88789798e+1*1e4,
                peak_wavelength=1.44888598e+0,
                peak_spectral_radiance =1.31061590e+1*1e4,
                band_radiance=2.88510571e+1*1e4,
                spectrum=read_temperature_data("2000.csv") )                                                  

)

"""
    testing_data(tag::String)

    returns benchmark data for specified tag
"""
function testing_data(tag::String)
    return benchmark_data[tag]
end
"""
    temperatures()

Returns numerical array of all temperatures in benchmark
"""
function temperatures()
    return [b.temperature for b in values(benchmark_data)]
end

end
