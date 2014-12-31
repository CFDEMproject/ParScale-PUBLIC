# Install/unInstall package files for ParScale
# mode = 0/1/2 for uninstall/install/update

mode=$1
# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# all package files with no dependencies

for file in *.cpp *.h; do
  action $file
done

# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then
  echo "...installing a new coupling model to style_coupling_model.h"
  if (test -e ../style_coupling_model.h) then
    sed -i '/liggghts/d' ../style_coupling_model.h  #remove old text related to ligghts
    echo '#include "coupling_model_liggghts.h"' >> ../style_coupling_model.h
  fi

elif (test $1 = 0) then
  if (test -e ../style_coupling_model.h) then
      echo "...cleaning style_coupling_model.h"
      sed -i '/liggghts/d' ../style_coupling_model.h  #remove old text related to ligghts
  fi

#  if (test -e ../Makefile.package.settings) then
#    sed -i -e '/^include.*pascal.*$/d' ../Makefile.package.settings
#  fi

fi
