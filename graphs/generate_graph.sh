cat graph.txt | dot -Tsvg -Kfdp -o viz_fdp.svg
cat graph.txt | dot -Tsvg -Kcirco -o viz_circo.svg
cat graph.txt | dot -Tsvg -o viz_dot.svg