#ifndef _DT_H_INCLUDED
#define _DT_H_INCLUDED

typedef struct DGNodeInfo_s {
  int id = -1;
  int feature_len = 0;
  int inDegree = 0;
  int outDegree = 0;
  int inArc[SMALL_BLOCK_SIZE] = {0};
  int outArc[SMALL_BLOCK_SIZE] = {0};
  char name[SMALL_BLOCK_SIZE];

} DGNodeInfo;

template <typename GraphData>
class Graph
{
public:
  Graph(
    string const name_, size_t numNodes_, GraphData &data_):
    m_name(name_),
    m_nodes(numNodes_, dash::BLOCKED),
    m_data(data_)
  {}

  dash::Array<DGNodeInfo> & nodes() {
    return m_nodes;
  }

  GraphData & data() {
    return m_data;
  }

  string const name() const {
    return m_name;
  }
private:
  string const m_name;
  dash::Array<DGNodeInfo> m_nodes;
  GraphData &m_data;
};

#endif
