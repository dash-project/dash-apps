#ifndef BLOCK_PREFETCHER__H
#define BLOCK_PREFETCHER__H

#include <map>
#include "MatrixBlock.h"
#include "common.h"
#include "ExtraeInstrumentation.h"


template<typename MatrixT>
class BlockPrefetcher {

public:
  using value_t = typename MatrixT::value_type;
  using block_t = MatrixBlock<MatrixT>;
  using block_pos_t = std::pair<size_t, size_t>;


  BlockPrefetcher(MatrixT& matrix, size_t block_size, size_t num_blocks)
  : block_size(block_size), num_blocks(num_blocks), matrix(&matrix)
  {
    buffer = new value_t[block_size*block_size*num_blocks];
  }

  ~BlockPrefetcher() {
    if (buffer) {
      delete[] buffer;
      buffer = nullptr;
    }
  }

  value_t *
  get(size_t i, size_t j) {
    value_t *result;
    block_pos_t block_pos = std::make_pair(i, j);
    auto it = prefetch_blocks.find(block_pos);
    if (it == prefetch_blocks.end()) {
      block_t block(*matrix, i, j);
      if (block.is_local()) {
        prefetch_blocks.insert(
          std::make_pair(block_pos, block.lbegin()));
        result = block.lbegin();
      } else {

        auto block_pre = &buffer[next_buf_pos()];
        dash::tasks::async(
          [=]() mutable {
            EXTRAE_ENTER(EVENT_PREFETCH);
            block.fetch_async(block_pre);
            while (!block.test()) {
              EXTRAE_EXIT(EVENT_PREFETCH);
              dash::tasks::yield(-1);
              EXTRAE_ENTER(EVENT_PREFETCH);
            }
            EXTRAE_EXIT(EVENT_PREFETCH);
          },
          //DART_PRIO_HIGH,
          dash::tasks::in(block),
          dash::tasks::out(block_pre)
        );
        prefetch_blocks.insert(
          std::make_pair(block_pos, block_pre));
        result = block_pre;
      }
    } else {
      result = it->second;
    }
    return result;
  }

  void
  clear() {
    prefetch_blocks.clear();
  }

private:

  size_t next_buf_pos() {
    size_t res = pos * block_size * block_size;
    pos = (pos + 1) % num_blocks;
    return res;
  }


private:
  value_t *buffer = nullptr;
  std::map<block_pos_t, value_t*> prefetch_blocks;
  size_t pos = 0;
  size_t block_size = 0;
  size_t num_blocks;
  MatrixT *matrix;
};

#endif // BLOCK_PREFETCHER__H
