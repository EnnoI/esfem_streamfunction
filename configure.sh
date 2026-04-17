#!/usr/bin/env bash

VERSION="2.10"
BRANCH="releases/2.10"
DUNE_DIR=${HOME}/DUNE_DIR/${VERSION}
DUNE_CONTROL="${DUNE_DIR}/dune-common/bin/dunecontrol"


# create build-scripts
cat <<EOF > build_debug.sh
#!/usr/bin/env bash

DUNE_DIR="${DUNE_DIR}"
DUNE_CONTROL="\${DUNE_DIR}/dune-common/bin/dunecontrol"

DUNE_CONTROL_PATH=\${DUNE_DIR}:. \${DUNE_CONTROL} --opts=\${DUNE_DIR}/debug.opts --only=esfem_stream all
EOF
chmod +x build_debug.sh

cat <<EOF > build_release.sh
#!/usr/bin/env bash

DUNE_DIR="${DUNE_DIR}"
DUNE_CONTROL="\${DUNE_DIR}/dune-common/bin/dunecontrol"

DUNE_CONTROL_PATH=\${DUNE_DIR}:. \${DUNE_CONTROL} --opts=\${DUNE_DIR}/release.opts --only=esfem_stream all
EOF
chmod +x build_release.sh


mkdir -p ${DUNE_DIR}
cd ${DUNE_DIR}

# download all needed dune modules into folder DUNE_DIR
git clone --branch ${BRANCH} https://gitlab.dune-project.org/core/dune-common.git
git clone --branch ${BRANCH} https://gitlab.dune-project.org/core/dune-geometry.git
git clone --branch ${BRANCH} https://gitlab.dune-project.org/core/dune-grid.git
git clone --branch ${BRANCH} https://gitlab.dune-project.org/core/dune-istl.git
git clone --branch ${BRANCH} https://gitlab.dune-project.org/core/dune-localfunctions.git
git clone --branch ${BRANCH} https://gitlab.dune-project.org/staging/dune-functions.git
git clone --branch ${BRANCH} https://gitlab.dune-project.org/staging/dune-typetree.git
git clone --branch ${BRANCH} https://gitlab.dune-project.org/extensions/dune-foamgrid.git
git clone --branch ${BRANCH} https://gitlab.dune-project.org/extensions/dune-alugrid.git
git clone --branch ${BRANCH} https://gitlab.dune-project.org/extensions/dune-curvedgeometry.git
git clone --branch ${BRANCH} https://gitlab.dune-project.org/extensions/dune-curvedgrid.git
git clone --branch master https://gitlab.dune-project.org/extensions/dune-vtk.git
git clone --branch ${BRANCH} https://gitlab.dune-project.org/extensions/dune-gmsh4.git
git clone --branch ${BRANCH} https://gitlab.mn.tu-dresden.de/mapo162b/amdis-surfaceextension.git
git clone --branch ${BRANCH} https://gitlab.dune-project.org/staging/dune-uggrid.git
git clone --branch ${BRANCH} https://gitlab.com/amdis/amdis.git

# create an opts file
cat <<EOF > debug.opts
CMAKE_FLAGS=" \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS=\"-DAMDIS_INFO_LEVEL=3 -DDUNE_ISTL_WITH_CHECKING=1 -DDUNE_CHECK_BOUNDS=1 -DDUNE_FMatrix_WITH_CHECKING=1 -DCHECK_RESERVEDVECTOR=1 -Wall\" \
  -DCMAKE_CXX_FLAGS_DEBUG=\"-O0 -g3 -fno-omit-frame-pointer\" \
  -DCMAKE_CXX_STANDARD=20 \
  -DDUNE_ENABLE_PYTHONBINDINGS:BOOL=0 \
  -DCMAKE_DISABLE_FIND_PACKAGE_MPI:BOOL=0 \
  -DCMAKE_DISABLE_FIND_PACKAGE_fmt=OFF \
  -DAMDIS_USE_SYSTEM_FMT=ON"
MAKE_FLAGS="-j8"
BUILDDIR="build-debug"
EOF

cat <<EOF > release.opts
CMAKE_FLAGS=" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS=\"-march=native\" \
  -DCMAKE_CXX_STANDARD=20 \
  -DDUNE_ENABLE_PYTHONBINDINGS:BOOL=0 \
  -DCMAKE_DISABLE_FIND_PACKAGE_MPI:BOOL=0 \
  -DCMAKE_DISABLE_FIND_PACKAGE_fmt=OFF \
  -DAMDIS_USE_SYSTEM_FMT=ON"
MAKE_FLAGS="-j8"
BUILDDIR="build-release"
EOF

# configure and build all dune modules
${DUNE_CONTROL} --opts=debug.opts all
${DUNE_CONTROL} --opts=release.opts all
