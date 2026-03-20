

    def replay(self, reactions: list[tuple[int, str, list[list[int]]]], until_frame=-1):
        """
        Replays reactions
        = runs the modification algorithm based on reaction names and atom indices
        """
        for frame, rx_name, frags in reactions:
            if 0 <= until_frame <= frame:
                break
            # find rx
            rx = self.top.rx_by_name(rx_name)
            assert rx is not None, f"{rx} could not be found."
            # find frags
            frags = [
                self.top.frag_by_name_and_atoms(frag_name, frag)
                for frag_name, frag in zip(rx.__reactants, frags)
            ]
            assert all(frag is not None for frag in frags)
            # modifications
            self.top.modification(frame, [(frags, rx)])
        self.system.reinitialize()
        # note to self: we don't make constraints and __vsites whole again as reactions can't modify those
        # so if they are made whole when the sim is constructed that's enough

    def step(self, steps=1, xtc=True, dm=True):
        """
            Do a step of the following:
            - steps MD steps
            - D/M algorithm if dm is true
            - reinitialize system if reactions happened
            - write an XTC frame if xtc is True
            - display info to logs and screen, % info given by self.i and self.md_steps
            - update self.i
        """
        start_time = time.time()
        self.i += steps
        percent = self.i/self.md_steps*100 if self.md_steps > 0 else 100
        self.logger.info(f"step {self.i}")
        if steps > 0:
            self.logger.info(f"md_steps {steps}")
            self.logger.info("MD start")
            try:
                self.system.do_steps(steps)
            except mm.OpenMMException as e:
                # TODO insert this to all openmm calls that can do an exception in some elegant manner, e.g. in Context
                self.logger.error(f"!!! OpenMM Exception !!!\n{e}")
                raise e

            self.logger.info("MD finished")
        if self.md_steps > 0 and steps > 0:
            ns_so_far = self.dt_ns * self.i
            time_left = self.last_step_time * (self.md_steps - self.i)
            time_fmt = format_time(time_left)
            reporter_data = " ".join(filter(None, [r.interactive_line() for r in self.reporters]))
            sys.stdout.write(
                f"\033[2K\rstep {self.i}"
                f"({ns_so_far:.2f} ns, "
                f"{percent:.1f}%) "
                f"{time_fmt} {reporter_data}"
            )
        if dm:
            self.logger.info("Detection start")
            pos, box = self.system.get_positions()
            reactions = self.top.detection(self.i, box, pos)
            self.logger.info("Detection finished")
            if len(reactions) > 0:
                self.logger.info("Modification start")
                reactions = self.top.modification(self.i, reactions)
                self.reactions += len(reactions)
                self.logger.info("Modification finished")
            if len(reactions) > 0 or self.force_reinitialize:
                self.logger.info("Reinitialize start")
                # TODO fix redundant reinitializes with update parameter in context for softcore
                self.system.reinitialize(self.force_reinitialize)
                for rep in self.reporters:
                    rep.post_reinitialize(self.i)
                self.top.toggle_sc(reactions, True)
                self.system.reinitialize(self.force_reinitialize)
                for rep in self.reporters:
                    rep.post_sc_enable(self.i)
                if len(reactions) > 0 and self.daemon_integrator:
                    self.integrator.set_reactions(reactions, self.system, self.top, self.i)
                self.top.toggle_sc(reactions, False)
                self.system.reinitialize(self.force_reinitialize)
                for rep in self.reporters:
                    rep.post_sc_disable(self.i)
                self.logger.info("Reinitialize finished")
            self.logger.info(f"reactions {len(reactions)}")

        if xtc:
            self.logger.info("XTC write start")
            self.system.write_xtc_frame(self.i, self.xtc_freq)
            self.logger.info("XTC write finished")
        end_time = time.time()
        # timing info update
        if steps > 0:
            step_time = (end_time - start_time) / steps
            if self.last_step_time > 0.:
                self.last_step_time = step_time * 0.01 + self.last_step_time * 0.99
            elif not xtc or not dm:
                # step 0 tends to have both xtc and dm as True, and is
                # usually unrepresentatively slow
                scale = self.dm_freq / self.xtc_freq
                if scale > 1.:
                    scale = 1. / scale
                if self.first_step_time == 0.:
                    # continuations might not start with an expensive step
                    scale = 0.
                self.last_step_time = step_time * (1 - scale) + scale * self.first_step_time
            else:
                # step 0 probably
                self.first_step_time = step_time
