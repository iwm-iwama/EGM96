#!/usr/bin/env ruby
#coding:utf-8

$VERSION = "iwm20210711"
# <<History>>
#   iwm20210606/20210711 Support TSV only.
#   iwm20210520 Support Fixed length, TSV, CSV.
#   iwm20200206 Support TSV only.
#   iwm20040819 Support Fixed length, TSV, CSV.

#-------------------------------------------------------------------------------
# This Ruby script, intpt.rb, was written based on a FORTRAN program, INTPT.F,
#   and by Yoshiyuki Iwama(iwm-iwama), August, 2004 and June, 2021.
#
# The program, INTPT.F, was provided by Professor Richard H. Rapp of The Ohio
#   State University, December, 1996.
#   It was put into its present form by the National Imagery and Mapping Agency
#   (NIMA), December, 1996.
#       https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84
#
#-------------------------------------------------------------------------------
# Ruby Script 'intpt.rb'
#   WGS 84 EGM96 15-Minute Calculator
#
# <<Environment>>
#   This script executing needs 'Ruby'.
#       http://www.ruby-lang.org/
#
# <<Execute>>
#   > ruby intpt.rb
#
# <<Option example>>
#   > ruby intpt.rb -g="ww15mgh.grd.tsv" -i="input.tsv" -o="outintpt.tsv"
#
#-------------------------------------------------------------------------------
# <<File Information>>
#   Input and Output data formats support TSV(Tab Separated Value).
#
#   <<Example>>
#     1.Input data
#       Latitude	Longitude	Remarks
#         38.628155	269.779155
#        -14.621217	305.021114
#        -90.000000	360.000000
#         35.360555	138.727222	Mt.Fuji seismograph
#
#     2.Output data
#       Latitude	Longitude	Elevation(m)	Remarks
#         38.628155	269.779155	-31.628
#        -14.621217	305.021114	-2.969
#        -90.000000	360.000000	999999.000
#         35.360555	138.727222	41.300	Mt.Fuji seismograph
#
#-------------------------------------------------------------------------------

$IFN = {
	"Grid File"   => "ww15mgh.grd.tsv",
	"Input File"  => "input.tsv"
}

$OFN = {
	"Output File" => "outintpt.tsv"
}

$LN72 = "------------------------------------------------------------------------"

#-------------------------------------------------------------------------------

# Speed Up!!
GC.disable

Signal.trap(:INT) do
	exit
end

IWINDO    = 4
NBDR      = IWINDO * 2  # 8
NLON      = 1441 + NBDR # 1449
NLON_NBDR = NLON - NBDR # 1441
NLAT      = 721
DLAT      = 0.25
DLON      = 0.25

$South    = 0.0
$North    = 0.0
$West     = 0.0
$East     = 0.0
$Dphi     = 0.0
$Dlam     = 0.0

$AryR     = []
$HashH    = {}

def HashKeyMake(i1, i2)
	# [1..721]*10000 + [1..1449]
	return (i1.to_i * 10000) + i2.to_i
end

def INITSP(aryY)
	#-------------
	# Initialize
	#-------------
	$AryR = Array.new((IWINDO + 1), 0.0)
	aryQ  = Array.new((IWINDO + 1), 0.0)

	#------------------
	# [2..IWINDO - 1]
	#------------------
	i1 = 2
	while i1 <= IWINDO - 1
		i2 = aryQ[i1 - 1] / 2.0 + 2
		aryQ[i1] = -0.5 / i2
			y1 = aryY[i1 + 1].to_f
			y2 = aryY[i1].to_f
			y3 = aryY[i1 - 1].to_f
		$AryR[i1] = (3 * (y1 - 2 * y2 + y3) - $AryR[i1 - 1] / 2.0) / i2
		i1 += 1
	end

	#------------------
	# Recalculation
	# [2..IWINDO - 1]
	#------------------
	i1 = 2
	while i1 <= (IWINDO - 1)
		$AryR[i1] = aryQ[i1] * $AryR[i1 + 1] + $AryR[i1]
		i1 += 1
	end

	aryQ = []
end

def INTERP(dmin, phi, dla)
	aryY  = []
	aryHC = []

	rho    = 57.29577951
	rearth = 6371000.0

	f1 = dmin * 1000 * rho
	f2 = $West - NBDR / 2 * $Dlam

	ilim = f1 / (rearth * DLAT)
	jlim = f1 / (rearth * DLON * Math.cos(($South + DLAT * NLAT / 2.0) / rho))

	ri = (phi - $South) / DLAT
	rj = (dla - f2) / DLON

	i6 = ri.floor - IWINDO / 2 + 1
	i7 = rj.floor - IWINDO / 2 + 1
	i8 = i6 + IWINDO - 1
	i9 = i7 + IWINDO - 1

	if i6 < 0 || i8 >= NLAT || i7 < 0 || i9 >= NLON
		return 999999.0
	end

	if i6 < ilim || i8 > (NLAT - ilim) || i7 < jlim || i9 > (NLON - jlim)
		return 999999.0
	end

	i1 = 1
	while i1 <= IWINDO
		i2 = 1
		while i2 <= IWINDO
			s = HashKeyMake(i6 + i1, i7 + i2)
				aryY[i2] = $HashH[s]
			i2 += 1
		end
		INITSP(aryY)
		aryHC[i1] = SPLINE(rj - i7 + 1, aryY)
		i1 += 1
	end
	INITSP(aryHC)
	valint = SPLINE(ri - i6 + 1, aryHC)

	aryY = []

	return valint
end

def BILIN(ri, rj)
	i1 = ri.floor
	i2 = rj.floor

	rn = ri - i1
	re = rj - i2

	if i1 < 1
		i1 = 1
		rn = 0.0
	elsif i1 >= NLAT
		i1 = NLAT - 1
		rn = 1.0
	end

	if i2 < 1
		i2 = 1
		re = 0.0
	elsif i2 >= NLON
		i2 = NLON - 1
		re = 1.0
	end

	rnm1 = 1 - rn
	rem1 = 1 - re

	s1 = HashKeyMake(i1, i2)
	s2 = HashKeyMake(i1 + 1, i2)
	s3 = HashKeyMake(i1, i2 + 1)
	s4 = HashKeyMake(i1 + 1, i2 + 1)

	return (rnm1 * rem1 * $HashH[s1]) + (rn * rem1 * $HashH[s2]) + (rnm1 * re * $HashH[s3]) + (rn * re * $HashH[s4])
end

def SPLINE(x, aryY)
	rtn = 0

	if x < 1
		rtn = aryY[1] + ((x - 1) * (aryY[2] - aryY[1] - $AryR[2] / 6.0))
	elsif x > IWINDO
		rtn = aryY[IWINDO] + (x - IWINDO) * (aryY[IWINDO] - aryY[IWINDO - 1] + $AryR[IWINDO - 1] / 6.0)
	else
		i1 = x.floor
		i2 = x - i1
			y1 = aryY[i1 + 1].to_f
			y2 = aryY[i1].to_f
			r1 = $AryR[i1 + 1].to_f
			r2 = $AryR[i1].to_f
		rtn = y2 + i2 * ((y1 - y2 - r2 / 3 - r1 / 6) + i2 * (r2 / 2 + i2 * (r1 - r2) / 6))
	end

	return rtn
end

def Line2Ary(ln)
	# TSV
	return ln.split("\t")

	## TSV | Fixed | CSV
	#	return ln.split(/\s+|,/)
end

def GrdF_read(grdFn)
	puts "\e[0;93mLoading a Grid File..."

	# Header [0..5]    # 6
	# Data   [0..1440] # 1441

	# lstrip(Remove leading whitespace)
	aryGrd = File.binread(grdFn).split("\n").map(&:lstrip)

	#---------
	# Header
	#---------
	$South, $North, $West, $East, $Dphi, $Dlam = Line2Ary(aryGrd[0])[0..5].map(&:to_f)

	#-----------------------------
	# Data Array 1441 * Line 721
	#-----------------------------
	j1 = NBDR / 2       # 4
	j2 = NLON_NBDR - j1 # 1437
	j3 = NLON - j1      # 1445

	$HashH = {}
	i2 = 0
	cnt = 1
	while cnt < aryGrd.size
		#
		# L181(Data=721) = L20 * 9 + L1
		#
		ary = Line2Ary(aryGrd[cnt])
		#
		# [5..1445]
		#
		i3 = 1
		while i3 <= NLON_NBDR
			s = HashKeyMake(NLAT - i2, i3 + j1)
				$HashH[s] = ary[i3 - 1]
			i3 += 1
		end
		#
		# [1..4]
		# [1446..1449]
		#
		i3 = 1
		while i3 <= j1
			s = HashKeyMake(NLAT - i2, i3)
				$HashH[s] = ary[j2 + i3 - 1]
			s = HashKeyMake(NLAT - i2, i3 + j3)
				$HashH[s] = ary[i3 - 1]
			i3 += 1
		end

		i2 += 1
		cnt += 1
	end
end

def main()
	ARGV.each do
		|s|
		if s =~ /^-g=/
			$IFN["Grid File"] = s.split("=")[1]
		elsif s =~ /^-i=/
			$IFN["Input File"] = s.split("=")[1]
		elsif s =~ /^-o=/
			$OFN["Output File"] = s.split("=")[1]
		end
	end

	# Clear screen
	system "cls || clear"

	print "\e[0;92m"
	puts $LN72
	printf(
		"> %s -g=\"%s\" -i=\"%s\" -o=\"%s\"\n",
		File.basename($0),
		$IFN["Grid File"],
		$IFN["Input File"],
		$OFN["Output File"]
	)
	puts $LN72
	print "\e[0;96m"

	iErr = 0

	# Input File, "ok" "NG"
	$IFN.each do
		|key, value|
		s = ""
		if File.exist?(value)
			s = "ok"
		else
			s = "\e[0;91mNG\e[0;96m"
			iErr = 1
		end
		printf("[%s] %-11s : '%s'\n", s, key, value)
	end

	# Output File, All "ok"
	$OFN.each do
		|key, value|
		printf("[%s] %-11s : '%s'\n", "ok", key, value)
	end

	if iErr > 0
		puts "\e[0;91mExit on error.\e[0;99m"
		exit
	end

	GrdF_read($IFN["Grid File"])

	rtn = ""
	outputCnt = 0

	File.open($IFN["Input File"], "rb") do
		|fs|
		fs.read.split("\n") do
			|ln|
			ln.strip!

			if ln.size > 0
				a1 = Line2Ary(ln)
				flat, flon = a1[0].to_f, a1[1].to_f
				m = INTERP(12.0, flat, flon)

				# TSV
				rtn << "#{flat}\t#{flon}\t#{m.floor(3)}"
				rtn << (a1[2] ? "\t#{a1[2]}\n" : "\n")

				## Fixed
				#	rtn << sprintf(
				#		"%14.7f%14.7f%12.3f",
				#		flat,
				#		flon,
				#		INTERP(12.0, flat, flon)
				#	)
				#	rtn << (a1[2] ? "  #{a1[2]}\n" : "\n")

				outputCnt += 1
			end
		end
	end

	print "\e[0;93m"
	puts "Output > #{outputCnt}"

	File.open($OFN["Output File"], "wb") do
		|fs|
		fs.write rtn
	end

	print "\e[0;99m"
end

main()
