#!/bin/bash
#SBATCH -t 05:00:00
#SBATCH -p fat
#SBATCH -N 1
cd /home/dijkbio/andre/CESM/src/run
module load pre2019
module load eb
module load Miniconda3
source activate CESM

# (
# sleep 0
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2000
# python run_OHC_integrals.py rcp 2001
# python run_OHC_integrals.py rcp 2002
# python run_OHC_integrals.py rcp 2003
# python run_OHC_integrals.py rcp 2004
# )&
# (
# sleep 20
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2005
# python run_OHC_integrals.py rcp 2006
# python run_OHC_integrals.py rcp 2007
# python run_OHC_integrals.py rcp 2008
# python run_OHC_integrals.py rcp 2009
# )&
# (
# sleep 40
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2010
# python run_OHC_integrals.py rcp 2011
# python run_OHC_integrals.py rcp 2012
# python run_OHC_integrals.py rcp 2013
# python run_OHC_integrals.py rcp 2014
# )&
# (
# sleep 60
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2015
# python run_OHC_integrals.py rcp 2016
# python run_OHC_integrals.py rcp 2017
# python run_OHC_integrals.py rcp 2018
# python run_OHC_integrals.py rcp 2019
# )&
# (
# sleep 80
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2020
# python run_OHC_integrals.py rcp 2021
# python run_OHC_integrals.py rcp 2022
# python run_OHC_integrals.py rcp 2023
# python run_OHC_integrals.py rcp 2024
# )&
# (
# sleep 100
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2025
# python run_OHC_integrals.py rcp 2026
# python run_OHC_integrals.py rcp 2027
# python run_OHC_integrals.py rcp 2028
# python run_OHC_integrals.py rcp 2029
# )&
# (
# sleep 120
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2030
# python run_OHC_integrals.py rcp 2031
# python run_OHC_integrals.py rcp 2032
# python run_OHC_integrals.py rcp 2033
# python run_OHC_integrals.py rcp 2034
# )&
# (
# sleep 140
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2035
# python run_OHC_integrals.py rcp 2036
# python run_OHC_integrals.py rcp 2037
# python run_OHC_integrals.py rcp 2038
# python run_OHC_integrals.py rcp 2039
# )&
# (
# sleep 160
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2040
# python run_OHC_integrals.py rcp 2041
# python run_OHC_integrals.py rcp 2042
# python run_OHC_integrals.py rcp 2043
# python run_OHC_integrals.py rcp 2044
# )&
# (
# sleep 180
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2045
# python run_OHC_integrals.py rcp 2046
# python run_OHC_integrals.py rcp 2047
# python run_OHC_integrals.py rcp 2048
# python run_OHC_integrals.py rcp 2049
# )&
# (
# sleep 200
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2050
# python run_OHC_integrals.py rcp 2051
# python run_OHC_integrals.py rcp 2052
# python run_OHC_integrals.py rcp 2053
# python run_OHC_integrals.py rcp 2054
# )&
# (
# sleep 220
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2055
# python run_OHC_integrals.py rcp 2056
# python run_OHC_integrals.py rcp 2057
# python run_OHC_integrals.py rcp 2058
# python run_OHC_integrals.py rcp 2059
# )&
# (
# sleep 240
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2060
# python run_OHC_integrals.py rcp 2061
# python run_OHC_integrals.py rcp 2062
# python run_OHC_integrals.py rcp 2063
# python run_OHC_integrals.py rcp 2064
# )&
# (
# sleep 260
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2065
# python run_OHC_integrals.py rcp 2066
# python run_OHC_integrals.py rcp 2067
# python run_OHC_integrals.py rcp 2068
# python run_OHC_integrals.py rcp 2069
# )&
# (
# sleep 280
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2070
# python run_OHC_integrals.py rcp 2071
# python run_OHC_integrals.py rcp 2072
# python run_OHC_integrals.py rcp 2073
# python run_OHC_integrals.py rcp 2074
# )&
# (
# sleep 300
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2075
# python run_OHC_integrals.py rcp 2076
# python run_OHC_integrals.py rcp 2077
# python run_OHC_integrals.py rcp 2078
# python run_OHC_integrals.py rcp 2079
# )&
# (
# sleep 320
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2080
# python run_OHC_integrals.py rcp 2081
# python run_OHC_integrals.py rcp 2082
# python run_OHC_integrals.py rcp 2083
# python run_OHC_integrals.py rcp 2084
# )&
# (
# sleep 340
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2085
# python run_OHC_integrals.py rcp 2086
# python run_OHC_integrals.py rcp 2087
# python run_OHC_integrals.py rcp 2088
# python run_OHC_integrals.py rcp 2089
# )&
# (
# sleep 360
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2090
# python run_OHC_integrals.py rcp 2091
# python run_OHC_integrals.py rcp 2092
# python run_OHC_integrals.py rcp 2093
# python run_OHC_integrals.py rcp 2094
# )&
# (
# sleep 380
# now=$(date +"%T")
# echo "Current time : $now"
# python run_OHC_integrals.py rcp 2095
# python run_OHC_integrals.py rcp 2096
# python run_OHC_integrals.py rcp 2097
# python run_OHC_integrals.py rcp 2098
# python run_OHC_integrals.py rcp 2099
# )&
# wait

# for i in np.arange(17): # for rcp: 20, 5*i
#     print('(')
#     print(f'sleep {i*20}')
#     print(f'python run_OHC_integrals.py hq {2000+3*i}')
#     print(f'python run_OHC_integrals.py hq {2000+3*i+1}')
#     print(f'python run_OHC_integrals.py hq {2000+3*i+2}')
# #     print(f'python run_OHC_integrals.py hq {2000+3*i+3}')
#     print(')&')
# print('wait')

(
sleep 0
# python run_OHC_integrals.py hq 2000
python run_OHC_integrals.py hq 2001
python run_OHC_integrals.py hq 2002
)&
(
sleep 20
python run_OHC_integrals.py hq 2003
python run_OHC_integrals.py hq 2004
python run_OHC_integrals.py hq 2005
)&
(
sleep 40
python run_OHC_integrals.py hq 2006
python run_OHC_integrals.py hq 2007
python run_OHC_integrals.py hq 2008
)&
(
sleep 60
python run_OHC_integrals.py hq 2009
python run_OHC_integrals.py hq 2010
python run_OHC_integrals.py hq 2011
)&
(
sleep 80
python run_OHC_integrals.py hq 2012
python run_OHC_integrals.py hq 2013
python run_OHC_integrals.py hq 2014
)&
(
sleep 100
python run_OHC_integrals.py hq 2015
python run_OHC_integrals.py hq 2016
python run_OHC_integrals.py hq 2017
)&
(
sleep 120
python run_OHC_integrals.py hq 2018
python run_OHC_integrals.py hq 2019
python run_OHC_integrals.py hq 2020
)&
(
sleep 140
python run_OHC_integrals.py hq 2021
python run_OHC_integrals.py hq 2022
python run_OHC_integrals.py hq 2023
)&
(
sleep 160
python run_OHC_integrals.py hq 2024
python run_OHC_integrals.py hq 2025
python run_OHC_integrals.py hq 2026
)&
(
sleep 180
python run_OHC_integrals.py hq 2027
python run_OHC_integrals.py hq 2028
python run_OHC_integrals.py hq 2029
)&
(
sleep 200
python run_OHC_integrals.py hq 2030
python run_OHC_integrals.py hq 2031
python run_OHC_integrals.py hq 2032
)&
(
sleep 220
python run_OHC_integrals.py hq 2033
python run_OHC_integrals.py hq 2034
python run_OHC_integrals.py hq 2035
)&
(
sleep 240
python run_OHC_integrals.py hq 2036
python run_OHC_integrals.py hq 2037
python run_OHC_integrals.py hq 2038
)&
(
sleep 260
python run_OHC_integrals.py hq 2039
python run_OHC_integrals.py hq 2040
python run_OHC_integrals.py hq 2041
)&
(
sleep 280
python run_OHC_integrals.py hq 2042
python run_OHC_integrals.py hq 2043
python run_OHC_integrals.py hq 2044
)&
(
sleep 300
python run_OHC_integrals.py hq 2045
python run_OHC_integrals.py hq 2046
python run_OHC_integrals.py hq 2047
)&
(
sleep 320
python run_OHC_integrals.py hq 2048
python run_OHC_integrals.py hq 2049
python run_OHC_integrals.py hq 2050
)&
wait


# python run_OHC_integrals.py ctrl  # 45 mins
# python run_OHC_integrals.py lpd   # 26 mins