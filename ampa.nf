#!/bin/env nextflow

/*
 * Copyright (c) 2013, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'AMPA-NF'.
 *
 *   Piper-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Piper-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with AMPA-NF.  If not, see <http://www.gnu.org/licenses/>.
 */


/* 
 * Authors: 
 * - Irantzu Anzar <iranmdl15@gmail.com>
 * - Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

params.in = 'example.fa'
params.out = 'bigampa.txt'


fastaFile = file(params.in)
if( !fastaFile.exists() ) {
  exit 1, "The specified input file does not exist: ${params.in}"
}


seq = channel()
fastaFile.chunkFasta { seq << it }

ampaOut = channel()
task ('ampa') {

    //  defines the Input and Output
    input '-':seq
    output '-':ampaOut

    // The BASH script to be executed - for each - sequence
    """
    cat - > input.file && AMPA-BIGTABLE.pl -in=input.file -noplot -rf=result -df=data
    cat input.file | grep '>' > /dev/stdout
    cat result | grep '#' > /dev/stdout
    """

}

def result = file(params.out)
println "Saving result at: ${result}"

ampaOut.eachWithIndex { str, index ->

  def lines = str.trim().split('\n')
  if( lines.size() != 2 ) {
      println "ERROR > Invalid AMPA ($index) result:\n${str}\n\n"
      return
  }

  def id = getIDs(lines[0])
  def val = getValues(lines[1])

  result << "${id[0]}\t${id[1]}\t${val[0]}\t${val[1]}\t${val[2]}\t${val[3]}\n"
}


def getIDs( line ) {

  def matcher = line =~ />(\S+).+gene:(\S+).*/
  if( matcher.matches() ) {
    def seqId = matcher[0][1]
    def geneId = matcher[0][2]
    return [seqId, geneId]
  }
  else {
    println "Bad ID line: $line"
    return []
  }

}


/*
 *  return the values in the following order 
 *  - stretches 
 *  - protLen: 
 *  - ampLen: 
 *  - propensity
 */
def getValues(result) {

 def rm = result =~ /# This protein has (\d+) bactericidal stretches and it has (\d+) amino acids. AMP length: (\d+) Best AMP Propensity: ([0-9\.]+)/

  if( rm.matches() ) {
    return [rm[0][1], rm[0][2], rm[0][3], rm[0][4]]
  }
  else {
    println "Bad result line: $result"
    return []   
  }
}


