.PHONY: all, clean

all: problem_statement.pdf

problem_statement.pdf: problem_statement.tex
	latexmk -pdf problem_statement

clean: 
	rm problem_statement-1_0.pdf problem_statement-1.asy problem_statement-1.pre problem_statement-1.tex problem_statement.aux problem_statement.fdb_latexmk problem_statement.fls problem_statement.log problem_statement.pdf problem_statement.pre

%.html : %.md
	Rscript -e "templater::render_template(md.file='$<', output='$@')"
