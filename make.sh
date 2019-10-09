#!/usr/bin/env bash
Rscript build.R _source/tutorial.Rmd;
Rscript build.R _source/advanced_tutorial.Rmd;
rm -rf docs
jekyll build
rm -rf docs
cp -r _site docs
rm docs/setup.md
jekyll serve
