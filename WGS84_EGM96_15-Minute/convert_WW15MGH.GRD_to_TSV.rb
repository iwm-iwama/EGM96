#!/usr/bin/env ruby
#coding:binary

$VERSION = "iwm20210606"
# <<History>>
#   iwm20210606
#   iwm20200206
#   iwm20200202

#-------------------------------------------------------------------------------
# Convert
#   'WW15MGH.GRD' of the WGS 84 EGM96 15-Minute Geoid Height File,
#       https://earth-info.nga.mil/php/download.php?file=egm-96interpolation
# To
#   'ww15mgh.grd.tsv' of TSV format
#
# <<Environment>>
#   This script executing needs 'Ruby'.
#       http://www.ruby-lang.org/
#
#-------------------------------------------------------------------------------

$IFN = "./WW15MGH.GRD"
$OFN = "./ww15mgh.grd.tsv"

#-------------------------------------------------------------------------------

Signal.trap(:INT) do
	exit
end

puts

#--------
# Input
#--------
$flg = 0

Dir.glob(File.dirname($IFN) + "/*") do
	|s|
	if s.downcase == $IFN.downcase
		$IFN = s
		$flg = 1
		break
	end
end

if $flg == 0
	printf("Not exist '%s'\n", $IFN)
	puts
	exit
end

#----------
# Convert
#----------
$rtn = ""

printf("%-4s '%s'\n", "From", $IFN)
File.open($IFN, "rb") do
	|fs|
	#
	# Header
	#
	while ln = fs.gets
		ln.strip!
		if ln.size > 0
			$rtn << ln.gsub(/\s+/, "\t") << "\n"
			break
		end
	end
	#
	# Data
	#
	cnt = 0
	while ln = fs.gets
		ln.strip!
		if ln.size > 0
			$rtn << ln.gsub(/\s+/, "\t")
			if (cnt += 1) < 181
				$rtn << "\t"
			else
				$rtn << "\n"
				cnt = 0
			end
		end
	end
end

printf("%-4s '%s'\n", "To", $OFN)
File.open($OFN, "wb") do
	|fs|
	fs.write $rtn
end

puts
