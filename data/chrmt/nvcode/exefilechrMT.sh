for filename in block1chrMT.sh block2chrMT.sh block3chrMT.sh block4chrMT.sh block5chrMT.sh block6chrMT.sh block7chrMT.sh block8chrMT.sh block9chrMT.sh block10chrMT.sh block11chrMT.sh block12chrMT.sh block13chrMT.sh block14chrMT.sh block15chrMT.sh block16chrMT.sh block17chrMT.sh block18chrMT.sh block19chrMT.sh block20chrMT.sh block21chrMT.sh 
do
 nohup ./$filename > $filename.log 2>&1 &
done
