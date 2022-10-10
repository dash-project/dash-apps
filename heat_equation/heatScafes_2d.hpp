#include <array>

#include "ScaFES.hpp"

constexpr size_t INIT_WIDTH = 5;

template<typename CT, std::size_t DIM>
class HeatEqnFDM : public ScaFES::Problem<HeatEqnFDM<CT,DIM>, CT, DIM> {
  public:
    /** All fields which are related to the underlying problem
     * are added in terms of an entry of the parameters of
     * type \c std::vector.
     * @param params Set of ScaFES parameters.
     * @param gg Global grid.
     * @param useLeapfrog Should the leap frog scheme be used?
     * @param nameDatafield Name of the fields.
     * @param stencilWidth Stencil width of the fields.
     * @param isKnownDf Is the data field are known or unknown one?
     * @param nLayers Number of layers at the global boundary.
     * @param defaultValue Default value of fields.
     * @param writeToFile How often should the data field be written to file.
     * @param computeError Should the Linf error between the numerical
     *                     and exact solution be computed?
     * @param geomparamsInit Initial guess of geometrical parameters.
     */
    HeatEqnFDM(ScaFES::Parameters const& params,
               ScaFES::GridGlobal<DIM> const& gg,
               bool useLeapfrog,
               std::vector<std::string> const& nameDatafield,
               std::vector<int> const& stencilWidth,
               std::vector<bool> const& isKnownDf,
               std::vector<int> const& nLayers = std::vector<int>(),
               std::vector<CT> const& defaultValue = std::vector<CT>(),
               std::vector<ScaFES::WriteHowOften> const& writeToFile
                 = std::vector<ScaFES::WriteHowOften>(),
               std::vector<bool> const& computeError = std::vector<bool>(),
               std::vector<CT> const& geomparamsInit = std::vector<CT>() )
        : ScaFES::Problem<HeatEqnFDM<CT, DIM>, CT, DIM>(params, gg, useLeapfrog,
                                                        nameDatafield, stencilWidth,
                                                        isKnownDf, nLayers,
                                                        defaultValue, writeToFile,
                                                        computeError, geomparamsInit)
        {
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        }

    /** Evaluates all fields at one given global inner grid node.
     *  @param vNew Set of all fields.
     *  @param idxNode Index of given grid node.
     */
    void evalInner(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& /*timestep*/) {
    }
    /** Evaluates all fields at one given global border grid node.
     *  @param vNew Set of all fields.
     *  @param idxNode Index of given grid node.
     *  @param timestep Given time step.
     */
    void evalBorder(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                    ScaFES::Ntuple<int,DIM> const& idxNode,
                    int const& timestep) {
    }

    /** Initializes all unknown fields at one given global inner grid node.
     *  @param vNew Set of all unknown fields (return value).
     *  @param idxNode Index of given grid node.
     */
    template<typename TT>
    void initInner(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                   std::vector<TT> const& /*vOld*/,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& timestep) {
        const auto& partition = this->globalGrid().partition(rank);
        auto first = partition.idxNodeFirstSub();
        if((idxNode[0] -  first[0]) < INIT_WIDTH && (idxNode[1] -  first[1]) < INIT_WIDTH) {
            vNew[0](idxNode) = 1;
        } else {
            vNew[0](idxNode) = 0;
        }
    }
    /** Initializes all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields (return value).
     *  @param vOld Set of all given fields.
     *  @param idxNode Index of given grid node.
     *  @param timestep Given tim step.
     */
    template<typename TT>
    void initBorder(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                    std::vector<TT> const& vOld,
                    ScaFES::Ntuple<int,DIM> const& idxNode,
                    int const& timestep) {
        this->template initInner<TT>(vNew, vOld, idxNode, timestep);
    }

    /** Updates all unknown fields at one given global inner grid node.
     *  @param vNew Set of all unknown fields at new time step (return value).
     *  @param vOld Set of all unknown fields at old time step.
     *  @param idxNode Index of given grid node.
     */
    template<typename TT>
    void updateInner(std::vector<ScaFES::DataField<TT,DIM>>& vNew,
                     std::vector<ScaFES::DataField<TT,DIM>> const& vOld,
                     ScaFES::Ntuple<int,DIM> const& idxNode,
                     int const& timestep) {
      auto center = vOld[0](idxNode);
      double dtheta = (vOld[0](this->connect(idxNode, 0)) + vOld[0](this->connect(idxNode, 1)) - 2 * center) / (dx * dx) +
                      (vOld[0](this->connect(idxNode, 2)) + vOld[0](this->connect(idxNode, 3)) - 2 * center) / (dy * dy);

        vNew[0](idxNode) = center + k * dtheta * dt;
    }
    /** Updates all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields at new time step (return value).
     *  @param idxNode Index of given grid node.
     */
    template<typename TT>
    void updateBorder(std::vector<ScaFES::DataField<TT,DIM>>& vNew,
                      std::vector<ScaFES::DataField<TT,DIM>>const& vOld,
                      ScaFES::Ntuple<int,DIM> const& idxNode,
                      int const& timestep) {
        auto center = vOld[0](idxNode);
        std::array<TT,4> stencil_values{};
        for(size_t i = 0; i < stencil_values.size(); ++i) {
            auto idxNode_tmp = this->connect(idxNode, i);
            if(idxNode_tmp[0] < 0 || idxNode_tmp < 0) {
                continue;
            }
            stencil_values[i] = vOld[0](idxNode_tmp);
        }
        double dtheta = (stencil_values[0] + stencil_values[1] - 2 * center) / (dx * dx) +
                        (stencil_values[2] + stencil_values[3] - 2 * center) / (dy * dy);

        vNew[0](idxNode) = center + k * dtheta * dt;
    }

    /** Updates (2nd cycle) all unknown fields at one given global inner grid node.
     *  \remarks Only important if leap frog scheme is used.
     */
    template<typename TT>
    void updateInner2(std::vector<ScaFES::DataField<TT,DIM>>&,
                      std::vector<ScaFES::DataField<TT,DIM>> const&,
                      ScaFES::Ntuple<int,DIM> const&,
                      int const&) { }

    /** Updates (2nd cycle) all unknown fields at one given global border
     *  grid node.
     *  \remarks Only important if leap frog scheme is used.
     */
    template<typename TT>
    void updateBorder2(std::vector<ScaFES::DataField<TT,DIM>>&,
                       std::vector<ScaFES::DataField<TT,DIM>>const&,
                       ScaFES::Ntuple<int,DIM> const&,
                       int const&) { }

    template<typename TT, size_t DIMENSION>
    double calc_energy(const std::vector<ScaFES::DataField<TT,DIMENSION>>& v) {
        const auto& partition = this->globalGrid().partition(rank);
        auto first = partition.idxNodeFirstSub();
        auto last = partition.idxNodeLastSub();
        TT energy{};
        for(int i = first[0]; i <= last[0]; ++i) {
            for(int j = first[1]; j <= last[1]; ++j) {
                energy += v[0](ScaFES::Ntuple<int,DIM>(i,j));
            }
        }

        double energy_final = 0;
        MPI_Reduce(&energy, &energy_final, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        return energy_final;
    }
private:
    static constexpr double dx = 1.0;
    static constexpr double dy = 1.0;
    static constexpr double dt = 0.05;
    static constexpr double k = 1.0;

    int rank = 0;
};
