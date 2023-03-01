import { pluginRegistry } from 'visyn_core/plugin';
import reg from './phovea';

/**
 * build a registry by registering all phovea modules
 */
// other modules
import 'visyn_core/phovea_registry';

// self
pluginRegistry.register('reaction-cime', reg);
