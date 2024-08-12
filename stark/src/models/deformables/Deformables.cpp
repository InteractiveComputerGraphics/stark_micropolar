#include "Deformables.h"

#include "../../utils/include.h"

stark::Deformables::Deformables(core::Stark& stark, spPointDynamics dyn)
	: point_sets(dyn)
{
	this->output = std::make_shared<DeformablesMeshOutput>(stark, dyn);
	this->lumped_inertia = std::make_shared<EnergyLumpedInertia>(stark, dyn);
	this->prescribed_positions = std::make_shared<EnergyPrescribedPositions>(stark, dyn);
	this->segment_strain = std::make_shared<EnergySegmentStrain>(stark, dyn);
	this->triangle_strain = stark.settings.models.enable_default_tri_strain ?  std::make_shared<EnergyTriangleStrain>(stark, dyn) : nullptr;
	this->discrete_shells = std::make_shared<EnergyDiscreteShells>(stark, dyn);
	this->tet_strain = stark.settings.models.enable_default_tet_strain ? std::make_shared<EnergyTetStrain>(stark, dyn) : nullptr;

    this->strain_kim_20 = std::make_shared<EnergyTriangleStrainKim20>(stark, dyn);
    this->strain_wen_23 = stark.settings.models.enable_model_wen23 ? std::make_shared<EnergyTriangleStrainWen23>(stark, dyn) : nullptr;
    this->strain_micropolar_shells = stark.settings.models.enable_model_mp_shell ? std::make_shared<EnergyMicropolarShells>(stark, dyn) : nullptr;
}

