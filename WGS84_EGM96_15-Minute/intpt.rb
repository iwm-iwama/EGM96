#!ruby
#coding:utf-8

$VERSION = ["iwm20200202", "iwm20040819"]
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
# This Ruby script, intpt.rb, was written based on a FORTRAN program, INTPT.F,
#   and by Yoshiyuki Iwama(iwm-iwama), August, 2004 and Janualy, 2020.
#
# The program, INTPT.F, was provided by Professor Richard H. Rapp of The Ohio
#   State University, December, 1996.
#   It was put into its present form by the National Imagery and Mapping Agency
#   (NIMA), December, 1996.
#       https://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
#
#-------------------------------------------------------------------------------
# <<File Information>>
#
#   1.FORTRAN Program INTPT.F
#       https://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/intpt.f
#
#   2.Geoid Height File
#       https://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/ww15mgh.grd.z
#
#   3.Test Input Data
#       https://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/input.dat
#
#       <<Example>> TSV format
#               38.6281550	269.7791550
#               -14.6212170	305.0211140
#               -90.0000000	360.0000000
#
#   4.Test Output Data
#       <<Example>> TSV format
#               38.6281550	269.7791550	-31.628
#               -14.6212170	305.0211140	-2.969
#               -90.0000000	360.0000000	999999.000	Err
#
#-------------------------------------------------------------------------------

def INTERP(dmin, phi, dla)
	aryY  = []
	aryHC = []

	rho    = 57.29577951
	rearth = 6371000.0

	f1 = dmin * 1000 * rho
	f2 = $west - NBDR / 2 * $dlam

	ilim = f1 / (rearth * $dlat)
	jlim = f1 / (rearth * $dlon * Math.cos(($south + $dlat * NLAT / 2.0) / rho))

	ri = (phi - $south) / $dlat
	rj = (dla - f2) / $dlon

	i6 = IFRAC(ri) - IWINDO / 2 + 1
	i7 = IFRAC(rj) - IWINDO / 2 + 1
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
	i1 = IFRAC(ri)
	i2 = IFRAC(rj)

	rn = ri - i1
	re = rj - i2

	case i1
		when i1 < 1
			i1 = 1
			rn = 0.0
		when i1 >= NLAT
			i1 = NLAT - 1
			rn = 1.0
	end

	case i2
		when i2 < 1
			i2 = 1
			re = 0.0
		when i2 >= NLON
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

def INITSP(aryY)
	#-------------
	# Initialize
	#-------------
	$AryR = [0.0] * (IWINDO + 1)
	aryQ  = [0.0] * (IWINDO + 1)

	#--------------------
	# [2 .. IWINDO - 1]
	#--------------------
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

	#--------------------
	# Recalculation
	# [2 .. IWINDO - 1]
	#--------------------
	i1 = IWINDO - 1
	while i1 >= 2
		$AryR[i1] = aryQ[i1] * $AryR[i1 + 1] + $AryR[i1]
		i1 -= 1
	end

	aryQ = []
end

def SPLINE(x, aryY)
	rtn = 0
	case x
		when x < 1
			rtn = aryY[1] + (
				(x - 1) * (aryY[2] - aryY[1] - $AryR[2] / 6.0)
			)
		when x > IWINDO
			rtn = aryY[IWINDO] + (x - IWINDO) * (aryY[IWINDO] - aryY[IWINDO - 1] + $AryR[IWINDO - 1] / 6.0)
		else
			i1 = IFRAC(x)
			i2 = x - i1
				y1 = aryY[i1 + 1].to_f
				y2 = aryY[i1].to_f
				r1 = $AryR[i1 + 1].to_f
				r2 = $AryR[i1].to_f
			rtn = y2 + i2 * (
				(y1 - y2 - r2 / 3 - r1 / 6) + i2 * (r2 / 2 + i2 * (r1 - r2) / 6)
			)
	end
	return rtn
end

def IFRAC(r)
	return r.floor
end

def HashKeyMake(i1, i2)
	return i1.to_s << "-" << i2.to_s
end

def Line2Ary(ln)
	# TSV only, fast.
	return ln.split("\t")

	# Compatible with the 'WW15MGH.GRD' file, but, slow about 20 percent.
	## return ln.split(/\s+/)
end

def GrdF_read(aGrdFn)
	File.open(aGrdFn, "rb") do
		|fs|

		# Header [0 .. 5]
		# Data   [0 .. 1440]
		ary = []

		#---------
		# Header
		#---------
		ln = fs.gets
			ln.strip!
		ary = Line2Ary(ln).map do
			|s|
			s.to_f
		end
			$south, $north, $west, $east, $dphi, $dlam = ary
		ary = []

		#-----------------------------
		# Data Array 1441 * Line 721
		#-----------------------------
		puts "                   << Loading a Grid File >>                    "
		puts "0%                            50%                           100%"
		puts "+-----------------------------+-----------------------------+-  "
		print "*"

		iGrh = 0

		j1 = NBDR / 2       # 4
		j2 = NLON_NBDR - j1 # 1437
		j3 = NLON - j1      # 1445

		$HashH = {}
		str = ""
		i1 = 0
		i2 = 0
		while ln = fs.gets
			ln.strip!
			if ln.size > 0
				str << ln << "\t"
				#
				# L181(Data=721) = L20 * 9 + L1
				#
				if (i1 += 1) >= 181
					i1 = 0
					ary = Line2Ary(str)
					str = ""
					#
					# [5 .. 1445]
					#
					i3 = 1
					while i3 <= NLON_NBDR
						s = HashKeyMake(NLAT - i2, i3 + j1)
							$HashH[s] = ary[i3 - 1]
						i3 += 1
					end
					#
					# [1 .. 4]
					# [1446 .. 1449]
					#
					i3 = 1
					while i3 <= j1
						s = HashKeyMake(NLAT - i2, i3)
							$HashH[s] = ary[j2 + i3 - 1]
						s = HashKeyMake(NLAT - i2, i3 + j3)
							$HashH[s] = ary[i3 - 1]
						i3 += 1
					end

					if (iGrh += 1) == 12
						iGrh = 0
						print "=>\b"
					end

					ary = []

					i2 += 1
				end
			end
		end
	end
	print "\n\n"
end

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------
BgnTime = Time.new

$IFN = {
	"Grid File"   => "ww15mgh.grd.tsv",
	"Input File"  => "input.tsv"
}

$OFN = {
	"Output File" => "outintpt.tsv"
}

ARGV.each do
	|s|
	case s
		when /^-g=/
			$IFN["Grid File"] = s.split("=")[1]

		when /^-i=/
			$IFN["Input File"] = s.split("=")[1]

		when /^-o=/
			$OFN["Output File"] = s.split("=")[1]
	end
end

66.times{print "-"}
puts
printf(
	"> %s -g=\"%s\" -i=\"%s\" -o=\"%s\"\n",
	File.basename($0),
	$IFN["Grid File"],
	$IFN["Input File"],
	$OFN["Output File"]
)
66.times{print "-"}
puts

iErr = 0

# Input File, "ok" "NG"
$IFN.each do
	|key, value|
	s = ""
	if File.exist?(value)
		s = "ok"
	else
		s = "NG"
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
	puts
	puts "Exit on error."
	puts
	exit
end
puts

#---------------
# Array / Hash
#---------------
$AryR = []
$HashH = {}

#------
# Int
#------
IWINDO    = 4
NBDR      = IWINDO * 2  # 8
NLON      = 1441 + NBDR # 1449
NLON_NBDR = NLON - NBDR # 1441
NLAT      = 721

#--------
# Float
#--------
$south = 0.0
$north = 0.0
$west  = 0.0
$east  = 0.0
$dphi  = 0.0
$dlam  = 0.0
$dlat  = 0.25
$dlon  = 0.25

#---------
# String
#---------
xmin = -50
xmax =  50
rdx  =   2

GrdF_read($IFN["Grid File"])

62.times{print "-"}
puts

rtn = ""

File.open($IFN["Input File"], "rb") do
	|fs|
	fs.each_line do
		|ln|
		ln.strip!
		if ln.size > 0
			ary = Line2Ary(ln).map do
				|s|
				s.to_f
			end
			flat, flon = ary

			un = INTERP(12.0, flat, flon)

			rtn << sprintf("%.7f\t%.7f\t%.3f", flat, flon, un)
			rtn << "\tErr" if un.abs == 999999.0
			rtn << "\n"
		end
	end
end

File.open($OFN["Output File"], "wb") do
	|fs|
	fs.write rtn
	print rtn
end

62.times{print "-"}
puts

printf("%.3fsec\n", Time.new - BgnTime)

puts
exit
