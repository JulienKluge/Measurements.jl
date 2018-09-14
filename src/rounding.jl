### rounding.jl
#
# Copyright (C) 2016, 2017 Mosè Giordano.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: uncertainty, error propagation, physics
#
# This file is a part of Measurements.jl.
#
# License is MIT "Expat".
#
### Commentary:
#
# This file contains rounding functions for the proper display of
# values with uncertainty within the measurement struct.
#
#
### Code:

export roundedmeasurementstr

export MeasurementDisplayMode
export DisplayTwoDigits, DisplayOneDigit, DisplayScientific, DisplayScientificConst, DisplayFull

#=
	flterrorresistant_roundup(num::T, sdigits::Int, flterroracceptance = 1.0e-14)
	roundedvalue

Returns an up-rounded value of `num` with `sdigits` significant digits. It ignores floating
point digits, often produced by rounding errors with a relative significance of `flterroracceptance`.
=#
function flterrorresistant_roundup(num::T, sdigits::Int, flterroracceptance = 1.0e-14) where {T <: AbstractFloat}
	rnum = round(num, RoundNearest, sigdigits = sdigits)
	relerror = (num - rnum) / (num + rnum)
	if abs(relerror) < flterroracceptance
		rnum
	else
		round(num, RoundUp, sigdigits = sdigits)
	end
end

#Display Mode type.
struct MeasurementDisplayMode{T} end

#Definition of all implemented types
const DisplayTwoDigits = MeasurementDisplayMode{:TwoDigits}()
const DisplayOneDigit = MeasurementDisplayMode{:OneDigit}()
const DisplayScientific = MeasurementDisplayMode{:Scientific}()
const DisplayScientificConst = MeasurementDisplayMode{:ScientificConst}()
const DisplayFull = MeasurementDisplayMode{:Full}()


#dirty hack to check if a number was formatted as exponential before padding
#TODO: find a way of forcing julia to format a float in a decimal string always
function rpad_notexp(str::String, count::Int) #rpad's with a string with "0", except it containing "e" to exclude exponentials
	if occursin("e", str)
		str
	else
		rpad(str, count, "0")
	end
end

#the power string to append
function composepowerstring(power::Int)
	power > 0  ? " * 10^$power" : " * 10^($power)"
end

#
# Formatting Functions for the respective display types
#

#(value ± error)
function composereg_rndstr(value::T, error::T, significantDgt::Int, firstValueDigit::Int, firstErrorDigit::Int, power::Int, negative::Bool) where T
	valueSignificantModifier = firstValueDigit - firstErrorDigit
	value = round(value, RoundNearestTiesUp, sigdigits = significantDgt + valueSignificantModifier)
	error = round(error, RoundUp, sigdigits = significantDgt) #errors are rounded up by default

	valueString = string(value)
	errorString = string(error)
	valueString = rpad_notexp(valueString, significantDgt + firstValueDigit + 2) #add trailing zeros
	errorString = rpad_notexp(errorString, significantDgt + firstErrorDigit + 2) #add trailing zeros

	signString = negative ? "-" : ""
	if power == 0
		"$signString($valueString ± $errorString)"
	else
		"$signString($valueString ± $errorString)$(composepowerstring(power))"
	end
end

#value(error)
function composeconst_rndstr(value::T, error::T, significantDgt::Int, firstValueDigit::Int, firstErrorDigit::Int, power::Int, negative::Bool) where T
	valueSignificantModifier = firstValueDigit - firstErrorDigit
	value = round(value, RoundNearestTiesUp, sigdigits = significantDgt + valueSignificantModifier)
	error = round(Int, error / 10.0^(firstErrorDigit - significantDgt + 1), RoundUp)

	valueString = string(value)
	errorString = string(error)
	valueString = rpad_notexp(valueString, significantDgt + abs(firstValueDigit) + 1) #add trailing zeros

	signString = negative ? "-" : ""
	if power == 0
		"$signString$valueString($errorString)"
	else
		"$signString$valueString($errorString)$(composepowerstring(power))"
	end
end

#(value ± error) not rounded
function composefull_rndstr(value::T, error::T, power::Int, negative::Bool) where T
	valueString = string(value)
	valueStrLen = length(valueString)
	errorString = string(error)
	errorStrLen = length(errorString)
	if valueStrLen > errorStrLen #add trailing zeros: error
		errorString = rpad_notexp(errorString, valueStrLen)
	elseif errorStrLen > valueStrLen  #add trailing zeros: value
		valueString = rpad_notexp(valueString, errorStrLen)
	end

	signString = negative ? "-" : ""
	if power == 0
		"$signString($valueString ± $errorString)"
	else
		"$signString($valueString ± $errorString)$(composepowerstring(power))"
	end
end
#
#
#

#main function
"""
	roundedmeasurementstring(value::T, error::T; exponentReducing::Int = 3, displaymode::MeasurementDisplayMode = DisplayTwoDigits)
	rounded_string_representation

Rounds and formats the value and error pair into the defined string representation, while optionally reducing the exponent.
"""
function roundedmeasurementstring(value::T, error::T; exponentReducing::Int = 3, displaymode::MeasurementDisplayMode = DisplayTwoDigits) where {T <: AbstractFloat}
	signFlag = value < zero(T) #get value sign and convert it to insertable string
	value = abs(value) #and remove sign from value

	firstValueDigit = value == 0 ? 0 : floor(Int, log10(value)) #get first nonzero digit from value
	firstErrorDigit = error == 0 ? 0 : floor(Int, log10(abs(error))) #get first nonzero digit from error
	#valueSignificantModifier = firstValueDigit - firstErrorDigit #expands or shrinks the significant value digits if the error lays in a different 10-log magnitude than the value
	
	powerDiv = 0
	if exponentReducing > 0
		powerDiv = div(firstValueDigit, exponentReducing) #magnitudes to store in exponents - reduce in 3-columns
		if powerDiv != 0 #reduce exponent
			powerDiv *= exponentReducing
			firstValueDigit -= powerDiv
			firstErrorDigit -= powerDiv
			powerDivModif = 10.0^powerDiv
			value /= powerDivModif
			error /= powerDivModif
		end
	end

	if displaymode == DisplayTwoDigits
		composedString = composereg_rndstr(value, error, 2, firstValueDigit, firstErrorDigit, powerDiv, signFlag)
	elseif displaymode == DisplayOneDigit
		composedString = composereg_rndstr(value, error, 1, firstValueDigit, firstErrorDigit, powerDiv, signFlag)
	elseif displaymode == DisplayScientific
		significantDgt = floor(Int, error / 10.0^firstErrorDigit) < 3 ? 2 : 1 #two dgt. if the error starts with 1 or 2 and one dgt. otherwise
		composedString = composereg_rndstr(value, error, significantDgt, firstValueDigit, firstErrorDigit, powerDiv, signFlag)
	elseif displaymode == DisplayScientificConst
		composedString = composeconst_rndstr(value, error, 2, firstValueDigit, firstErrorDigit, powerDiv, signFlag)
	else #MeasurementDisplayMode.DisplayFull
		composedString = composefull_rndstr(value, error, powerDiv, signFlag)
	end

	composedString
end

"""
	roundedmeasurementstring(m::Measurement{<:AbstractFloat}; exponentReducing::Int = 3, displaymode::MeasurementDisplayMode = DisplayTwoDigits)
	rounded_string_representation

Rounds and formats the value and error pair into the defined string representation, while optionally reducing the exponent.
"""
function roundedmeasurementstr(m::Measurement{<:AbstractFloat}; exponentReducing::Int = 3, displaymode::MeasurementDisplayMode = DisplayTwoDigits)
	roundedmeasurementstring(promote(m.val, m.err)..., exponentReducing = exponentReducing, displaymode = displaymode)
end