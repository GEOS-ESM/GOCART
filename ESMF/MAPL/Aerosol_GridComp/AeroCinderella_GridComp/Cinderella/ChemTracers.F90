module ChemTracers_mod
    use ChemTracers_private_mod, only: read_tracer_config, TracerMap

    implicit none
    private

    public read_tracer_config
    public TracerMap
end module ChemTracers_mod