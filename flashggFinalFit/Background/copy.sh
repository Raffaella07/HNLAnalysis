indirname="outdir_testMGRealData/bkgfTest-Data/"

# assumes eos mounted on t3home
target_dirname="testMGRealData"

rm -rf /eos/home-m/mratti/www/BHNL/recoAnalysis/fTest/$target_dirname
cp -r $indirname /eos/home-m/mratti/www/BHNL/recoAnalysis/fTest/$target_dirname
cp HTACCESS /eos/home-m/mratti/www/BHNL/recoAnalysis/fTest/$target_dirname/.htaccess

