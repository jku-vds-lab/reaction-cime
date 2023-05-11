import { createAsyncThunk, createEntityAdapter, createSlice, PayloadAction } from '@reduxjs/toolkit';
import { ISecureItem } from 'visyn_core/security';
import { ReactionCIMEBackendFromEnv } from '../Backend/ReactionCIMEBackend';

export type Project = { name: string; id: string; file_status: string } & ISecureItem;

const projectAdapter = createEntityAdapter<Project>();

const initialState = { projects: projectAdapter.getInitialState() };

export const deleteProject = createAsyncThunk('projects/delete', async (id: string, { dispatch }) => {
  await ReactionCIMEBackendFromEnv.deleteFile(id);
  return id;
});

export const projects = createSlice({
  name: 'projects',
  initialState,
  reducers: {
    addProject(state, action: PayloadAction<Project>) {
      projectAdapter.addOne(state.projects, action.payload);
    },
    setProjects(state, action: PayloadAction<Project[]>) {
      projectAdapter.setAll(state.projects, action.payload);
    },
  },
  extraReducers: (builder) => {
    builder.addCase(deleteProject.fulfilled, (state, action) => {
      projectAdapter.removeOne(state.projects, action.meta.arg);
    });
  },
});

export const ProjectActions = { ...projects.actions };

export const syncProjects = createAsyncThunk('projects/sync', async (_, { dispatch }) => {
  const files = await ReactionCIMEBackendFromEnv.getUploadedFiles();
  dispatch(projects.actions.setProjects(files));
});
