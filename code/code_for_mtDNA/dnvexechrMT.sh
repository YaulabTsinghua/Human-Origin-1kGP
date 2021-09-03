for filename in dnvblock1chrMT.sh dnvblock2chrMT.sh dnvblock3chrMT.sh dnvblock4chrMT.sh dnvblock5chrMT.sh dnvblock6chrMT.sh dnvblock7chrMT.sh dnvblock8chrMT.sh dnvblock9chrMT.sh dnvblock10chrMT.sh dnvblock11chrMT.sh dnvblock12chrMT.sh dnvblock13chrMT.sh dnvblock14chrMT.sh dnvblock15chrMT.sh dnvblock16chrMT.sh dnvblock17chrMT.sh dnvblock18chrMT.sh dnvblock19chrMT.sh dnvblock20chrMT.sh dnvblock21chrMT.sh 
do
 nohup ./$filename > $filename.log 2>&1 &
done
