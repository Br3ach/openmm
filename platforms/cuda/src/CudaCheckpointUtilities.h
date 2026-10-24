#ifndef OPENMM_CUDACHECKPOINTUTILITIES_H_
#define OPENMM_CUDACHECKPOINTUTILITIES_H_
/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012-2022 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */


#include "CudaContext.h"
#include "CudaArray.h"
#include "openmm/common/ContextSelector.h"
#include <algorithm>
#include <istream>
#include <limits>
#include <ostream>
#include <vector>

namespace OpenMM { namespace detail {
static const size_t CudaCheckpointChunkBytes = 8*1024*1024;
inline void checkCheckpointCuda(CUresult result, const CudaContext& context,
                                const ArrayInterface& array, size_t offset) {
    if (result != CUDA_SUCCESS)
        throw OpenMMException("CUDA checkpoint transfer for "+array.getName()+" at byte "+std::to_string(offset)+": "+context.getErrorString(result));
}
inline size_t checkpointArrayBytes(const ArrayInterface& array) {
    const int elementSize = array.getElementSize();
    if (elementSize <= 0 || array.getSize() > std::numeric_limits<size_t>::max()/elementSize)
        throw OpenMMException("Invalid CUDA checkpoint array size for "+array.getName());
    return array.getSize()*elementSize;
}
inline void transferCudaCheckpointArray(ArrayInterface& array, std::ostream* output, std::istream* input) {
    CudaContext& context = dynamic_cast<CudaContext&>(array.getContext());
    ContextSelector selector(context);
    const size_t bytes = checkpointArrayBytes(array);
    const CUdeviceptr pointer = context.unwrap(array).getDevicePointer();
    if (bytes > std::numeric_limits<CUdeviceptr>::max()-pointer)
        throw OpenMMException("CUDA checkpoint pointer overflow for "+array.getName());
    // A synchronous raw copy must not race work queued on a non-default stream.
    checkCheckpointCuda(cuStreamSynchronize(context.getCurrentStream()), context, array, 0);
    std::vector<char> buffer(std::min(bytes, CudaCheckpointChunkBytes));
    for (size_t offset = 0; offset < bytes; ) {
        const size_t count = std::min(bytes-offset, buffer.size());
        if (output != NULL) {
            checkCheckpointCuda(cuMemcpyDtoH(buffer.data(), pointer+offset, count), context, array, offset);
            output->write(buffer.data(), count);
            if (!*output)
                throw OpenMMException("Error writing CUDA checkpoint array "+array.getName());
        }
        else {
            input->read(buffer.data(), count);
            if (!*input)
                throw OpenMMException("Error reading CUDA checkpoint array "+array.getName());
            checkCheckpointCuda(cuMemcpyHtoD(pointer+offset, buffer.data(), count), context, array, offset);
        }
        offset += count;
    }
}
inline void writeCudaCheckpointArray(ArrayInterface& array, std::ostream& stream) {
    transferCudaCheckpointArray(array, &stream, NULL);
}
inline void readCudaCheckpointArray(ArrayInterface& array, std::istream& stream) {
    transferCudaCheckpointArray(array, NULL, &stream);
}
} }
#endif
