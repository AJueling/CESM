import sys
sys.path.append("..")
from bd_analysis_indices import AnalyzeIndex as AI
run = str(sys.argv[1])
extent = str(sys.argv[2])
AI().Pacific_EOF_analysis(run=run, extent=extent)