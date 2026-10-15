# Converting the floes to FixedWidthFloes, copying them to the device and back, and writing
# them back into the StructArray. timestep_sim! does this around the ported steps. It's
# the cost that goes away once the floes stay on the device between timesteps. The baseline
# has nothing like it.

function transfer!(floes, backend)
    dev_floes = adapt(backend, Subzero.FixedWidthFloes(floes))
    KernelAbstractions.synchronize(backend)
    host_floes = adapt(Array, dev_floes)
    Subzero.update_floes!(floes, host_floes)
    free!(dev_floes)
    return
end

function to_device!(fwf, backend)
    dev = adapt(backend, fwf)
    KernelAbstractions.synchronize(backend)
    free!(dev)
end

register_step!(;
    title = "host ↔ device: `FixedWidthFloes`, `adapt`, `update_floes!`",
    baseline = ctx -> nothing,
    current = function (ctx, backend, FT)
        floes = load_floes(ctx.n)
        fwf = Subzero.FixedWidthFloes(floes)
        # Every sample starts from the floes in fwf. This is the same as resetting them.
        trial = @bench(transfer!($floes, $backend), Subzero.update_floes!($floes, $fwf))
        record!(ctx, "round trip", trial; primary = true)
        trial = @bench(Subzero.FixedWidthFloes($floes))
        record!(ctx, "FixedWidthFloes(floes)", trial)
        trial = @bench(to_device!($fwf, $backend))
        record!(ctx, "adapt to device", trial)
        dev = adapt(backend, fwf)
        trial = @bench(adapt(Array, $dev))
        record!(ctx, "adapt to host", trial)
        free!(dev)
        trial = @bench(Subzero.update_floes!($floes, $fwf))
        record!(ctx, "update_floes!", trial)
    end,
)
